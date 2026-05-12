# Procedimento de Teste — Validação Modular do Framework eVTOL Tailsitter

> **Objetivo**: validar cada componente isoladamente ANTES de integrar o sistema completo. Cada teste tem critério de PASS/FAIL objetivo, tolerâncias declaradas e ação corretiva sugerida em caso de falha.

## Filosofia de teste

A revisão sênior identificou que a falha de tracking observada (~2000 m de erro) muito provavelmente NÃO é um bug isolado — pode ser uma combinação de:
- Mapeamento `force_to_attitude` mal-condicionado em hover (suspeita primária)
- Ganhos SO(3) sub-dimensionados para inércia de 1500 kg (suspeita secundária)
- Possíveis problemas dimensionais em alocador/cascata

A única forma de isolar a causa-raiz é **testar cada bloco em separado**, com sinais de entrada controlados, e só então integrar. Isto é o padrão DO-178C / MIL-HDBK-516 para sistemas críticos de voo.

## Como executar

Em MATLAB, dentro de `evtol_simulation/`:

```matlab
% Roda toda a suíte de testes
cd evtol_simulation
run_all_tests

% Ou um teste específico
cd tests
test_quat_utils
test_force_to_attitude   % o que mais nos interessa
```

Saída esperada: relatório formatado `[PASS] / [FAIL]` por teste, com sumário final.

## Estrutura da bateria

| Etapa | Suíte | Testes | Foco |
|---|---|---|---|
| 1a | `test_quat_utils` | 12 | Hamilton mul, conj, kinematic, exp/log |
| 1b | `test_so3_utils` | 8 | hat/vee, errMat, exp/log |
| 1c | `test_aircraft_dynamics` | 7 | 6-DOF, queda livre, conservação energia |
| 2a | `test_aerodynamics` | 8 | Viterna 360°, strip theory, slipstream |
| 2b | `test_propulsion` | 8 | BEMT, ESC dynamics, induced velocity |
| 3a | `test_so3_controller` | 6 | Estabilidade, ganho × inércia |
| 3b | `test_indi_controller` | 5 | Filtros, incremento |
| 3c | `test_nmpc` | 5 | Hover, step, restrição de tilt |
| 3d | `test_allocator` | 6 | Pseudoinversa, saturação, autoridade |
| 3e | `test_force_to_attitude` | 7 | **Mapeamento crítico (suspeita primária)** |
| 4a | `test_flatness` | 6 | Min-snap, BCs, derivadas |
| 4b | `test_trajectory_opt` | 3 | Custo energético, GA |
| 5a | `test_rk4` | 4 | Convergência, conservação |
| 5b | `test_dryden` | 4 | Variância, reset, estatística |
| 5c | `test_atmosphere` | 3 | ISA, monotonicidade |
| **Total** | **15 suítes** | **92 testes** | — |

---

## Etapa 1 — Cinemática e Dinâmica

### 1a. `test_quat_utils.m`

| # | Teste | Critério | Tolerância |
|---|---|---|---|
| Q1 | Identidade Hamilton: $q\otimes[1,0,0,0]=q$ | $\|q-q'\|<\epsilon$ | $1e-12$ |
| Q2 | Conjugado: $q\otimes q^* = [1,0,0,0]$ | $\|q\otimes q^*-e\|<\epsilon$ | $1e-12$ |
| Q3 | $R(I) = I_3$ | $\|R-I\|<\epsilon$ | $1e-12$ |
| Q4 | $\|R(q)\|=1$ (ortogonalidade) | $\|R^\top R-I\|<\epsilon$ | $1e-10$ |
| Q5 | Roundtrip Euler↔quat para $(\phi,\theta,\psi)\in[-\pi,\pi]^3$ random | erro angular < $\epsilon$ | $1e-9$ |
| Q6 | Cinemática: $\omega=0\Rightarrow\dot q=0$ | $\|\dot q\|<\epsilon$ | $1e-15$ |
| Q7 | Cinemática: $\omega=[0,0,1]$, $q=I\Rightarrow\dot q=\frac{1}{2}[0,0,0,1]$ | element-wise | $1e-12$ |
| Q8 | ExpMap(0) = identidade | $\|q-e\|<\epsilon$ | $1e-15$ |
| Q9 | LogMap(ExpMap($\phi$)) = $\phi$ para $\|\phi\|<1$ | $\|\phi-\phi'\|<\epsilon$ | $1e-9$ |
| Q10 | Rotação: $R_z(90°)\cdot[1,0,0]=[0,1,0]$ | $\|v-v'\|<\epsilon$ | $1e-12$ |
| Q11 | errMul: $q_e(q_d,q_d)=[1,0,0,0]$ | element-wise | $1e-12$ |
| Q12 | Norm preservation após múltiplas operações | $\|\|q\|-1\|<\epsilon$ | $1e-9$ |

**Diagnóstico**: falha em Q3-Q5 indica convenção de quatérnio errada (Hamilton vs JPL). Falha em Q6-Q7 quebra a integração de atitude — bug crítico.

### 1b. `test_so3_utils.m`

| # | Teste | Critério |
|---|---|---|
| S1 | $\text{vee}(\text{hat}(v))=v$ | erro $<1e-15$ |
| S2 | hat antissimétrica: $\text{hat}(v)+\text{hat}(v)^\top=0$ | $<1e-15$ |
| S3 | $\exp(0)=I$ | $<1e-15$ |
| S4 | $\exp(\log(R))=R$ para R aleatório | Frobenius $<1e-9$ |
| S5 | errMat($R,R$)=0 | $<1e-12$ |
| S6 | errMat: pequenos ângulos $\theta$ → $e_R\approx\theta\hat n$ | erro relativo $<1e-3$ |
| S7 | $\det(\exp(\phi))=1$ | $<1e-12$ |
| S8 | Rodrigues consistente com Euler-Rodrigues | $<1e-9$ |

### 1c. `test_aircraft_dynamics.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| D1 | Queda livre (sem F, sem M) | $\dot v_{NED}=g$ | $\|\dot v - [0;0;g]\| < 1e-12$ |
| D2 | Conservação norma quaternion após 1000 passos RK4 | hover, sem M | $\|\|q\|-1\| < 1e-6$ |
| D3 | Hover trim: F_prop = m·g (NED -Z) | M=0, V=0 | $\dot v_{NED}\approx 0$ |
| D4 | Pure rotation: M=[1;0;0] | F=0 | $\dot\omega = J^{-1}[1;0;0]$ |
| D5 | Acoplamento giroscópico: $\omega=[1;0;0]$, M=0 | $\dot p=v$ | trajetória esperada |
| D6 | $J_{xz}$ coupling: $\omega=[1;0;0]$ produz $\dot\omega_z\neq 0$ | sem M | sinal correto |
| D7 | Rotor gyro: 4 rotores Ω=200 rad/s, body $\omega=[0;1;0]$ | M_gyro $\neq 0$ | direção correta (precessão) |

**Diagnóstico crítico**: D1 falhando = bug de sinal de gravidade. D2 falhando = RK4 não está renormalizando q.

---

## Etapa 2 — Aero-Propulsão

### 2a. `test_aerodynamics.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| A1 | Sem ar parado: V=0 | strip a 0 m/s | F=M=0 |
| A2 | NACA-0015 a $\alpha=0$ | V=50 m/s | $C_L\approx 0$, $C_D=C_{D0}$ |
| A3 | Faixa linear: $C_L(\alpha)=C_{L_\alpha}\alpha$ | $\alpha\in[-10°,10°]$ | erro $<2\%$ |
| A4 | Estol: $C_L$ cai após $\alpha_s=15°$ | $\alpha=20°$ | $C_L(20°)<C_L(15°)$ |
| A5 | Placa plana 90°: $C_D\approx C_{D,\max}\approx 1.2$ | $\alpha=90°$ | $\|C_D-1.2\|<0.3$ |
| A6 | Simetria: $C_L(-\alpha)=-C_L(\alpha)$ | NACA-0015 simétrica | $\|C_L(\alpha)+C_L(-\alpha)\|<0.05$ |
| A7 | Pressão dinâmica $\propto V^2$ | dobrar V | F quadruplica |
| A8 | Slipstream eleva V_local | rotor empurrando, V_inf=0 | strip vê V > 0 |

### 2b. `test_propulsion.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| P1 | Ω=0 → T=0, Q=0 | comando zero | $T,Q<1e-9$ |
| P2 | Estática: $T=k_T\Omega^2$ | Ω = Ω_max/2 | erro $<1\%$ |
| P3 | Atuador 1ª ordem: convergência em $5\tau$ | $T_{cmd}=2000$ N degrau | $T_{actual}/T_{cmd}>0.99$ após $5\tau_m$ |
| P4 | Saturação: $T_{cmd}>T_{max}$ | comando $10^5$ N | $T_{actual}\le T_{max}$ |
| P5 | Velocidade induzida hover: $v_i=\sqrt{T/(2\rho A)}$ | T=3700 N, V_x=0 | erro $<2\%$ |
| P6 | Forward flight: $v_i$ diminui com $V_x$ | V_x = 30 m/s, mesmo T | $v_i(V=30)<v_i(0)$ |
| P7 | Reaction torque: $Q$ oposto ao spin | spin_dir = +1 | $M_z<0$ |
| P8 | 4 rotores em trim equalizam: $\sum T_i = m g$ | hover | erro $<5\%$ |

**Diagnóstico**: P3 falhando = $\tau_m$ errado, motor não acompanha. P5 falhando = momentum theory furada → slipstream errado.

---

## Etapa 3 — Controle

### 3a. `test_so3_controller.m` ⚠️ **CRÍTICO**

| # | Teste | Setup | Critério |
|---|---|---|---|
| C1 | Equilíbrio: $R=R_d$, $\omega=\omega_d=0$ | sem integral | $M\approx J\cdot 0 + \omega\times J\omega = 0$ |
| C2 | Erro pequeno: $R_d=R\cdot\exp(0.01\hat e_y)$ | $\omega=0$ | sinal de M correto |
| C3 | Bandwidth check: $\sqrt{k_R/J_{yy}}\ge 5$ rad/s | gains atuais | **PROVÁVEL FALHA** (gain × J = 220/3800 → ω_n=0.24 rad/s ≪ 5) |
| C4 | Damping ratio $\zeta\in[0.5,1.5]$ | atuais | verificar |
| C5 | Anti-windup: integral satura em I_max | step persistente | $\|e_I\|\le I_{\max}$ |
| C6 | Recuperação: erro 30°→0° em $<3\tau_n$ | $\omega_n=10$ rad/s nominal | tempo de assentamento |

**Critério de bandwidth (CRÍTICO)**:
$$\omega_n = \sqrt{k_R/J} \ge 5\text{ rad/s para tailsitter de massa moderada}$$

Com gains atuais e $J_{yy}=3800$: $\omega_n=\sqrt{220/3800}=0.24$ rad/s → **insuficiente, explica divergência**.

Gains corretos para $\omega_n=10$ rad/s:
$$k_R \approx J \cdot \omega_n^2,\quad k_\Omega \approx 2\zeta\omega_n J$$

### 3b. `test_indi_controller.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| N1 | Inicialização: filtros zerados | reset | estados $=0$ |
| N2 | Convergência LP: degrau em $\dot\omega$ | excitação degrau | $\dot\omega_{filt}$ converge em $1/\omega_{LP}$ |
| N3 | Estado estacionário: $q=q_d$, $\omega=0$ | sem dinâmica | $\Delta u\approx 0$ |
| N4 | $G^+$ bem-condicionado | $\sigma_{\min}/\sigma_{\max}>10^{-6}$ | número de cond. |
| N5 | Resposta a degrau de qd: produz $\Delta u\neq 0$ | qe não-trivial | direção correta |

### 3c. `test_nmpc.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| M1 | Hover hold: ref = current state, persistente | $V=0$ | $f_{cmd}\approx [0;0;-g]$ |
| M2 | Step longitudinal: ref +10 m em N | accel commanded > 0 | $f_{cmd,x}>0$ |
| M3 | Restrição de tilt obedecida | ref muito agressiva | $\arctan(f_{xy}/f_z)\le$ tilt_max |
| M4 | Saturação $a_{max}$ obedecida | distance >> | $\|f_{cmd}\|\le a_{max}$ |
| M5 | Warm-start: solução continua | dois passos consec. | redução de iterações |

### 3d. `test_allocator.m`

| # | Teste | Setup | Critério |
|---|---|---|---|
| L1 | Hover: $\nu=[0;0;-mg;0;0;0]$ | sem moments | $T_i\approx mg/4=3680$ N |
| L2 | Pitch moment puro | $\nu=[0;0;0;0;1000;0]$ | rotores L/R diferenciam |
| L3 | Yaw moment puro | $\nu=[0;0;0;0;0;500]$ | rotores spin assimétrico |
| L4 | Saturação top: $\nu$ excessivo | $T_{cmd}>T_{max}$ | $T_i$ saturado, residual $\neq 0$ |
| L5 | Rank do B (com cant): rank $\ge 4$ | numérico | autoridade em todos eixos |
| L6 | Allocator surface scaling: q_bar=0 surfaces 0 | hover | $\delta_e=0$ correto |

### 3e. `test_force_to_attitude.m` ⚠️ **CRÍTICO PRIMÁRIO**

| # | Teste | Setup | Critério Preciso |
|---|---|---|---|
| F1 | Hover puro: $F_{des}=[0;0;-mg]$, $\psi=0$ | tailsitter pointing up | $q_d$ corresponde a pitch=+90°, roll=0, yaw=0 |
| F2 | Cruise nivelado: $F_{des}=[0;0;-mg]+[Da;0;0]$ | accel longitudinal $a$ | $q_d$ tilta forward de $\arctan(a/g)$ |
| F3 | Continuidade: $\Delta F_{des}\to 0\Rightarrow\Delta q_d\to 0$ | finitas diferenças | Lipschitz |
| F4 | $\det(R(q_d))=+1$ | qualquer F | proper rotation |
| F5 | $\hat x_B(q_d)$ paralelo a $F_{des}$ | qualquer F | $\hat x_B \cdot \hat F = 1$ |
| F6 | Heading consistente: rotacionar $\psi_d$ rotaciona $\hat y_B$ no plano $\perp\hat x_B$ | $F_{des}=$const | $\hat y_B$ varia com $\psi_d$ |
| F7 | Sem singularidade: $\hat F\to[0;0;-1]$ | hover degenerado | sem NaN |

**Diagnóstico**: F1, F4 ou F5 falhando = ESTE é o bug primário. A função produz $q_d$ inconsistente, o loop interno persegue alvo errado, e a cascata diverge.

---

## Etapa 4 — Trajetória

### 4a. `test_flatness.m`

| # | Teste | Critério |
|---|---|---|
| T1 | $\sigma(0)=p_0$, $\sigma(T)=p_1$ | erro $<1e-9$ |
| T2 | $\dot\sigma(0)=\dot\sigma(T)=0$ (rest-to-rest) | erro $<1e-9$ |
| T3 | $\ddot\sigma(0)=\ddot\sigma(T)=0$ | erro $<1e-9$ |
| T4 | $\dddot\sigma(0)=\dddot\sigma(T)=0$ | erro $<1e-9$ |
| T5 | Continuidade $C^3$ entre segmentos | derivada $\le 3$ contínua |
| T6 | Recover: thrust direction aponta na diagonal correta de aceleração | $\hat x_B$ paralelo a $m(\ddot\sigma-g)$ | erro angular $<0.01°$ |

### 4b. `test_trajectory_opt.m`

| # | Teste | Critério |
|---|---|---|
| O1 | Baseline (Flat puro): tempo total = soma alocada | erro $<1e-12$ |
| O2 | GA reduz custo: $J_{ga}\le J_{baseline}$ | obrigatório se ativado |
| O3 | Sem violação de saturação no ótimo | $T_i\le T_{max}\forall t$ |

---

## Etapa 5 — Ambiente

### 5a. `test_rk4.m`

| # | Teste | Critério |
|---|---|---|
| R1 | Linear: $\dot x=Ax$ converge ordem 4 | erro $\propto h^4$ |
| R2 | Renorm de quaternion: $\|q\|\le 1+1e-9$ após 1000 passos | sempre |
| R3 | Conservação: oscilador harmônico ideal | erro de energia $<10^{-3}$ por período |
| R4 | Estabilidade A-region: pequeno passo respeita | $\|x\|$ limitado |

### 5b. `test_dryden.m`

| # | Teste | Critério |
|---|---|---|
| W1 | Disabled: retorna `constant_NED` | exato |
| W2 | Enabled: variância $\sigma^2$ próxima de teórica | $|\sigma_{emp}^2/\sigma_{theo}^2-1|<0.3$ (10000 amostras) |
| W3 | Reset: estados zerados | $x_u, x_v, x_w = 0$ |
| W4 | Saída zero-mean (após reset, longo run) | $\|\bar v\|<0.2\sigma$ |

### 5c. `test_atmosphere.m`

| # | Teste | Critério |
|---|---|---|
| Z1 | Sea level: $\rho(0)=1.225$ kg/m³ | erro $<0.001$ |
| Z2 | 11 km: $\rho<0.4$ kg/m³ | sim |
| Z3 | Monotonicidade: $\rho$ decrescente | sempre |

---

## Master Runner

`run_all_tests.m` executa toda a bateria, agrega contagens, e imprime relatório no formato:

```
==================================================
   eVTOL Tailsitter — Test Suite Report
==================================================
[STAGE 1] Cinemática/Dinâmica
  test_quat_utils ............ PASS (12/12)
  test_so3_utils ............. PASS (8/8)
  test_aircraft_dynamics ..... PASS (7/7)
[STAGE 2] Aero-Propulsão
  test_aerodynamics .......... PASS (8/8)
  test_propulsion ............ PASS (8/8)
[STAGE 3] Controle
  test_so3_controller ........ FAIL (4/6)  <-- C3, C6 (gain bandwidth)
  test_indi_controller ....... PASS (5/5)
  test_nmpc .................. PASS (5/5)
  test_allocator ............. PASS (6/6)
  test_force_to_attitude ..... FAIL (5/7)  <-- F1, F6 (hover singularity)
[STAGE 4] Trajetória
  test_flatness .............. PASS (6/6)
  test_trajectory_opt ........ PASS (3/3)
[STAGE 5] Ambiente
  test_rk4 ................... PASS (4/4)
  test_dryden ................ PASS (4/4)
  test_atmosphere ............ PASS (3/3)
==================================================
TOTAL: 84/92 tests passed (91.3%)
FAILURES: 8 tests in 2 suites
==================================================
```

## Critério de "Pronto para Integração"

O sistema só pode ser integrado (e a missão completa rodada) **DEPOIS** que todas as 15 suítes acusem PASS. Se alguma suíte falhar, a respectiva etapa precisa ser corrigida antes de prosseguir.

Em particular: **3a (SO(3))** e **3e (force_to_attitude)** são as suspeitas primárias para o erro de tracking observado. Esperamos que falhem na primeira execução, e a correção delas deve eliminar a divergência.

## Sequência recomendada de execução

1. Rodar `run_all_tests` → identificar falhas
2. Para cada suíte com FAIL, examinar o teste específico e a saída diagnóstica
3. Corrigir o módulo afetado
4. Re-rodar a suíte específica até PASS
5. Re-rodar `run_all_tests` para confirmar que correções não quebraram outros testes
6. Após 15/15 PASS: rodar `main.m` (missão completa) — esperar erro de tracking < 5 m em cruzeiro
