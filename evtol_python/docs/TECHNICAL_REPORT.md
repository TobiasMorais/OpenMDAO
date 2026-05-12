# Relatório Técnico Detalhado
## Modelagem e Simulação 6-DOF de eVTOL Tailsitter Biplace 1500 kg

**Versão**: 1.0 (Maio 2026)
**Plataforma**: MATLAB R2020a+ (OOP)
**Repositório**: github.com/TobiasMorais/OpenMDAO/tree/claude/evtol-simulation-framework-ARdun
**Total de linhas de código**: ~4500 (excluindo testes e docs)

---

## Sumário

1. [Resumo e Escopo](#1-resumo)
2. [Revisão de Literatura](#2-revisão-de-literatura)
3. [Descrição da Aeronave](#3-descrição-da-aeronave)
4. [Modelagem Cinemática e Dinâmica (6-DOF)](#4-modelagem-cinemática-e-dinâmica)
5. [Modelagem Aerodinâmica](#5-modelagem-aerodinâmica)
6. [Modelagem Propulsiva](#6-modelagem-propulsiva)
7. [Arquitetura de Controle](#7-arquitetura-de-controle)
8. [Planejamento de Trajetória](#8-planejamento-de-trajetória)
9. [Ambiente e Perturbações](#9-ambiente-e-perturbações)
10. [Métodos Numéricos](#10-métodos-numéricos)
11. [Tabela Completa de Parâmetros](#11-tabela-completa-de-parâmetros)
12. [Validação e Limites](#12-validação-e-limites)
13. [Referências](#13-referências)

---

## 1. Resumo

Este relatório documenta a modelagem física, matemática e implementação computacional de uma plataforma de simulação 6-DOF para veículos eVTOL não-convencionais do tipo **tailsitter biplace de 1500 kg com 4 rotores fixos**. A plataforma cobre cinemática rígida com quatérnios, dinâmica Newton-Euler com tensor de inércia assimétrico ($J_{xz}$), aerodinâmica Viterna-Corrigan 360°, modelo propulsivo BEMT com slipstream, atuadores de 1ª ordem, controle em cascata SO(3) + PD+FF, planejamento de trajetória por planicidade diferencial com refinamento genético, e ambiente atmosférico ISA com turbulência Dryden MIL-F-8785C.

**Faixa de validade do modelo**: Mach ≤ 0.3, altitudes ≤ 11 km (troposfera), ângulos de ataque ±90° (envelope completo via Viterna-Corrigan).

---

## 2. Revisão de Literatura

### 2.1 Dinâmica de aeronaves não-convencionais

A modelagem 6-DOF segue formulação clássica de **Stevens, Lewis & Johnson (2015)** [1] e **Etkin & Reid (1996)** [2], adaptada para aeronaves tipo tailsitter conforme **Lustosa et al. (2017)** [3]. O uso de quatérnios Hamilton segue **Markley & Crassidis (2014)** [4], evitando singularidades em pitch ±90° críticas para tailsitter.

### 2.2 Aerodinâmica em alto ângulo de ataque

O modelo de extrapolação pós-estol é o **modelo de Viterna & Corrigan (1981)** [5], originalmente desenvolvido para aerogeradores e extensivamente validado em aeronaves tilt-wing (**Anderson 1991** [6]). A teoria de tiras (strip theory) com slipstream-correção de propeller segue **Patterson & Bowman (2014)** [7] e **Borer et al. (2016)** [8] (X-57 Maxwell).

### 2.3 Modelagem propulsiva

A momentum theory clássica de Glauert (**Leishman 2006** [9]) é a base para velocidade induzida axial. O modelo de slipstream contraction (k_s = 1.0 conservador) segue **Drela (2014)** [10]. Atuadores ESC modelados como sistema de 1ª ordem (**Bristeau et al. 2009** [11]).

### 2.4 Controle em cascata

A malha interna SO(3) implementa a lei geométrica de **Lee, Leok & McClamroch (2010)** [12] — referência seminal para controle de quadrirrotores e UAVs vetorizados. A malha externa PD com feedforward de aceleração segue **Mellinger & Kumar (2011)** [13] (diferencial-flatness para quadrirrotores), comprovadamente robusta em rastreamento de trajetória.

Alternativa moderna INDI (Incremental Nonlinear Dynamic Inversion) baseada em **Smeur, Chu & de Croon (2016)** [14] está implementada e testada como malha interna selecionável.

### 2.5 Planejamento de trajetória

Polinômios mínimo-snap de 8ª ordem ($C^3$-contínuos) seguem **Mellinger, Mueller & Kumar (2012)** [15]. Refinamento por algoritmo genético implementa a otimização energética sugerida em **Bouabdallah, Murrieri & Siegwart (2004)** [16].

### 2.6 Turbulência atmosférica

O modelo de Dryden segue a especificação militar **MIL-F-8785C** (1980) [17], categoria de baixa altitude (<1000 ft). Atmosfera padrão ISA conforme **U.S. Standard Atmosphere (1976)** [18].

### 2.7 Aeronaves de referência

Os parâmetros físicos são derivados por escalamento de **Joby S4** [19] (massa 2200 kg, eVTOL urbano), **Beta Alia 250** [20] (eVTOL fixed-wing tilt-thrust) e dados de tailsitter de pesquisa **Verling et al. (2016)** [21].

---

## 3. Descrição da Aeronave

### 3.1 Configuração geral

| Item | Valor | Unidade |
|---|---|---|
| Tipo | Tailsitter biplace VTOL | — |
| Massa total (m) | 1500 | kg |
| Peso (W = mg) | 14 715 | N |
| Capacidade | 2 ocupantes | — |
| Massa da bateria | 420 | kg |
| Capacidade energética | 280 | kWh |

### 3.2 Geometria aerodinâmica

| Item | Valor | Unidade |
|---|---|---|
| Envergadura (b) | 11.0 | m |
| Corda aerodinâmica média (c) | 1.4 | m |
| Área da asa (S) | 15.0 | m² |
| Razão de aspecto (AR = b²/S) | 8.07 | — |
| Eficiência de Oswald (e) | 0.85 | — |
| Perfil aerodinâmico | NACA-0015 (simétrico para ops bidirecionais) | — |
| AoA de sustentação zero ($\alpha_0$) | 0 | rad |
| Inclinação de sustentação ($C_{L_\alpha}$) | 5.7 | 1/rad |
| $C_{L,max}$ | 1.35 | — |
| AoA de estol ($\alpha_{stall}$) | 15 | deg |
| Coef. arrasto parasítico ($C_{D_0}$) | 0.022 | — |
| $C_{m,0}$ | 0 | — |
| Derivada $C_{m_\alpha}$ | -0.45 | 1/rad |

### 3.3 Empenagens

| Estabilizador | Área (m²) | Braço (m) | $C_{L_\alpha}$ ou $C_{Y_\beta}$ |
|---|---|---|---|
| Horizontal | 2.5 | 4.5 | 4.5 1/rad |
| Vertical | 1.6 | 4.5 | -0.55 (CY_β) |

### 3.4 Superfícies de controle

| Superfície | Eficácia | Saturação |
|---|---|---|
| Elevator | 0.55 (dCL/dδe) | ±25° |
| Ailerons | 0.18 | ±25° |
| Rudder | 0.45 | ±25° |
| Servo τ | 40 ms | — |

### 3.5 Inércia

Tensor de inércia (plano de simetria X-Z, $J_{xy}=J_{yz}=0$):

$$J = \begin{bmatrix} J_{xx} & 0 & -J_{xz} \\ 0 & J_{yy} & 0 \\ -J_{xz} & 0 & J_{zz} \end{bmatrix} = \begin{bmatrix} 2200 & 0 & -180 \\ 0 & 3800 & 0 \\ -180 & 0 & 5500 \end{bmatrix} \text{ kg·m}^2$$

Centro de massa em coordenadas do corpo: $r_{CG} = [0.05; 0; 0]$ m.

### 3.6 Sistema propulsivo

**4 rotores fixos** dispostos aos pares em cada ponta de asa:

| Rotor | Posição (m) [x; y; z] | Pitch cant | Yaw cant | Spin |
|---|---|---|---|---|
| 1 (R, superior) | [0.30; +5.0; -0.15] | +2° | +1° | +1 (CCW) |
| 2 (R, inferior) | [0.30; +5.0; +0.15] | -2° | +1° | -1 (CW) |
| 3 (L, superior) | [0.30; -5.0; -0.15] | +2° | -1° | -1 (CW) |
| 4 (L, inferior) | [0.30; -5.0; +0.15] | -2° | -1° | +1 (CCW) |

**Parâmetros por rotor**:

| Item | Valor | Unidade |
|---|---|---|
| Raio (R) | 0.90 | m |
| Área do disco ($A = \pi R^2$) | 2.545 | m² |
| Número de pás | 5 | — |
| Corda da pá | 0.10 | m |
| Pitch coletivo das pás | 12° | — |
| Solidez ($\sigma = N_b c / \pi R$) | 0.177 | — |
| $C_{L_\alpha}$ da pá | 5.7 | 1/rad |
| $C_{D_0}$ da pá | 0.012 | — |
| Inércia rotor ($J_r$) | 0.32 | kg·m² |
| $\Omega_{max}$ | 230 (~2200 rpm) | rad/s |
| Empuxo máximo por rotor | 4500 | N |
| Coeficiente $k_T$ | $T_{max}/\Omega_{max}^2$ = 0.0850 | N·s² |
| Coeficiente $k_Q$ | $0.025 k_T R$ = 0.00191 | N·m·s² |
| Constante de tempo ESC ($\tau_m$) | 60 | ms |
| **Empuxo total disponível** | 18 000 | N |
| **T/W** | 1.22 | — |

---

## 4. Modelagem Cinemática e Dinâmica

### 4.1 Sistemas de referência

- **Frame Inercial NED** ($\mathcal{I}$): $X_I$=Norte, $Y_I$=Leste, $Z_I$=Down (paralelo a $\mathbf g$)
- **Frame Corpo** ($\mathcal{B}$): $X_B$=longitudinal, $Y_B$=ponta direita, $Z_B$=ventral

### 4.2 Vetor de estado (13 componentes)

$$\mathbf x = \begin{bmatrix} \mathbf p_{NED} \in \mathbb{R}^3 \\ \mathbf v_{NED} \in \mathbb{R}^3 \\ \mathbf q \in S^3 \\ \boldsymbol\omega^B \in \mathbb{R}^3 \end{bmatrix}$$

### 4.3 Quatérnios Hamilton

Convenção escalar-primeiro: $\mathbf q = [q_0, q_1, q_2, q_3]^\top$, $\|\mathbf q\| = 1$.

**Matriz de rotação corpo→mundo**:

$$R(\mathbf q) = (q_0^2 - \mathbf q_v^\top \mathbf q_v)\,I + 2\mathbf q_v \mathbf q_v^\top - 2 q_0\,[\mathbf q_v]_\times$$

onde $[\mathbf q_v]_\times$ é a matriz antissimétrica de $\mathbf q_v = [q_1, q_2, q_3]^\top$.

**Cinemática do quatérnio**:

$$\dot{\mathbf q} = \frac{1}{2}\Omega(\boldsymbol\omega^B)\,\mathbf q, \quad \Omega(\boldsymbol\omega) = \begin{bmatrix} 0 & -\omega_x & -\omega_y & -\omega_z \\ \omega_x & 0 & \omega_z & -\omega_y \\ \omega_y & -\omega_z & 0 & \omega_x \\ \omega_z & \omega_y & -\omega_x & 0 \end{bmatrix}$$

**Erro multiplicativo** (anti-unwinding):

$$\mathbf q_e = \mathbf q_d^{-1} \otimes \mathbf q, \quad \text{se } q_{e,0}<0 \Rightarrow \mathbf q_e \leftarrow -\mathbf q_e$$

### 4.4 Equações de Newton-Euler (corpo rígido)

**Translação** (integração no frame NED para facilidade de navegação):

$$\dot{\mathbf p}_{NED} = \mathbf v_{NED}$$
$$\dot{\mathbf v}_{NED} = \frac{1}{m}\,R(\mathbf q)\,(\mathbf F^B_{aero} + \mathbf F^B_{prop}) + \mathbf g_{NED}$$

com $\mathbf g_{NED} = [0; 0; 9.80665]^\top$ m/s².

**Rotação** (no frame do corpo, com $J_{xz}$ acoplado):

$$J\,\dot{\boldsymbol\omega}^B + \boldsymbol\omega^B \times (J\,\boldsymbol\omega^B) + \mathbf M^B_{gyro,rotor} = \mathbf M^B_{aero} + \mathbf M^B_{prop}$$

### 4.5 Torque giroscópico de rotores

$$\mathbf M^B_{gyro,rotor} = \boldsymbol\omega^B \times \sum_{i=1}^{N_r} J_r\,\Omega_i^{(s)}\,\hat{\mathbf e}_i$$

onde $\Omega_i^{(s)} = \mathrm{sign}_i \cdot \Omega_i$ (sinal por direção de rotação) e $\hat{\mathbf e}_i$ é o eixo de empuxo do rotor $i$ no corpo (com cant aplicado).

### 4.6 Álgebra SO(3) (controle interno)

**Operadores**:
- $\mathrm{hat}(\mathbf v) = [\mathbf v]_\times$ (skew-symmetric)
- $\mathrm{vee}(X) = \mathrm{hat}^{-1}(X)$ para $X$ antissimétrica
- $\exp_{SO(3)}(\boldsymbol\phi) = I + \sin(\theta)K + (1-\cos\theta)K^2$, $K = \boldsymbol\phi/\theta$ (Rodrigues)

**Erro de atitude geométrico** (Lee-Leok-McClamroch):

$$\mathbf e_R = \frac{1}{2}\mathrm{vee}(R_d^\top R - R^\top R_d)$$

### 4.7 Sanity checks numéricos

| Verificação | Resultado | Tolerância |
|---|---|---|
| Queda livre: $a_{NED} = [0;0;g]$ | erro $< 10^{-12}$ | $10^{-12}$ ✓ |
| Trim de hover $T = mg$ | $\|a\| < 10^{-9}$ | $10^{-9}$ ✓ |
| Conservação norma quaternion (1000 RK4) | erro $< 10^{-15}$ | $10^{-6}$ ✓ |
| Acoplamento $J_{xz}$: $\omega_x=1$ produz $\dot\omega_z \ne 0$ | sinal correto | — ✓ |

---

## 5. Modelagem Aerodinâmica

### 5.1 Decomposição em strips

A aeronave é decomposta em 4 superfícies portantes:
- Asa direita (`wing_R`): metade da área $S/2 = 7.5$ m², posição $[0.05; +b/4; 0]$
- Asa esquerda (`wing_L`): idêntica, $y$ invertido
- Estabilizador horizontal (`htail`): $S_t = 2.5$ m², posição $[-4.5; 0; 0]$
- Estabilizador vertical (`vtail`): $S_v = 1.6$ m², posição $[-4.5; 0; -0.5]$

Para cada strip $k$:

$$\mathbf V_{local,k} = \mathbf V^B + \boldsymbol\omega^B \times \mathbf r_k + \mathbf V_{slip,k}$$

### 5.2 Ângulo de ataque e sideslip

$$\alpha_k = \mathrm{atan2}(-V_{z,k}, V_{x,k})$$
$$\beta_k = \mathrm{asin}\left(\frac{V_{y,k}}{\|\mathbf V_{local,k}\|}\right)$$

### 5.3 Modelo Viterna-Corrigan (envelope 360°)

**Regime linear** ($|\alpha| \le \alpha_{stall}$):
$$C_L(\alpha) = C_{L_\alpha}(\alpha - \alpha_0)$$
$$C_D(\alpha) = C_{D_0} + \frac{C_L^2}{\pi AR \cdot e}$$

**Regime pós-estol** (extrapolação Viterna):

$$C_L(\alpha) = A_1 \sin(2\alpha) + A_2 \frac{\cos^2\alpha}{\sin\alpha}$$
$$C_D(\alpha) = B_1 \sin^2\alpha + B_2 \cos\alpha$$

onde:
$$C_{D,max} = 1.11 + 0.018\,AR = 1.255$$
$$A_1 = \frac{C_{D,max}}{2}, \quad A_2 = (C_{L,s} - C_{D,max}\sin\alpha_s\cos\alpha_s)\frac{\sin\alpha_s}{\cos^2\alpha_s}$$
$$B_1 = C_{D,max}, \quad B_2 = \frac{C_{D,s} - C_{D,max}\sin^2\alpha_s}{\cos\alpha_s}$$

Saturação dura final: $C_L \in [-2.0, +2.0]$, $C_D \in [0, 2.5]$ (defesa contra divisões near-zero).

### 5.4 Forças e momentos por strip

$$\bar q_k = \frac{1}{2}\rho \|\mathbf V_{local,k}\|^2$$
$$L_k = \bar q_k\,S_k\,C_L(\alpha_k), \quad D_k = \bar q_k\,S_k\,C_D(\alpha_k)$$

Para asas e htail (lift no plano X-Z do corpo):
$$\hat{\mathbf e}_{drag} = -\frac{\mathbf V_{local}}{\|\mathbf V_{local}\|}$$
$$\hat{\mathbf e}_{lift} = \frac{[-V_z; 0; V_x]}{\|[-V_z; 0; V_x]\|}$$

Para vtail (sideforce primário):
$$Y = \bar q\,S_v\,C_{Y_\beta}\,\beta$$

**Soma total**:
$$\mathbf F^B_{aero} = \sum_k \mathbf F_k, \quad \mathbf M^B_{aero} = \sum_k \mathbf r_k \times \mathbf F_k$$

### 5.5 Momento de arfagem estático

$$C_m = C_{m_0} + C_{m_\alpha}\,\alpha$$
$$M_{y,aero} \mathrel{+}= \bar q\,S\,c\,C_m$$

### 5.6 Contribuição das superfícies de controle

$$M_y \mathrel{+}= -\eta_{elev}\,S_{ht}\,L_{ht}\,\delta_e\,\bar q$$
$$M_x \mathrel{+}= \eta_{ail}\,S\,b\,\delta_a\,\bar q$$
$$M_z \mathrel{+}= -\eta_{rud}\,S_{vt}\,L_{vt}\,\delta_r\,\bar q$$

### 5.7 Ponto de máxima eficiência aerodinâmica

Para polar parabólica $C_D = C_{D_0} + K C_L^2$ com $K = 1/(\pi AR \cdot e) = 0.0467$:

$$C_{L,LDmax} = \sqrt{C_{D_0}/K} = 0.687$$
$$V_{LDmax} = \sqrt{\frac{2W}{\rho S \sqrt{C_{D_0}/K}}} = 48.2 \text{ m/s} = 174 \text{ km/h}$$
$$\alpha_{LDmax} = C_{L,LDmax}/C_{L_\alpha} = 6.9°$$
$$\left(\frac{L}{D}\right)_{max} = \frac{1}{2\sqrt{C_{D_0}K}} = 15.6$$

---

## 6. Modelagem Propulsiva

### 6.1 Eixo de empuxo com cant

Para cada rotor $i$ com cant pitch $\alpha_p^i$ e yaw $\alpha_y^i$:

$$\hat{\mathbf e}_i = R_z(\alpha_y^i)\,R_y(\alpha_p^i)\,\hat{\mathbf x}_B$$

### 6.2 Dinâmica do ESC (1ª ordem)

$$\tau_m\,\dot\Omega_i + \Omega_i = \Omega_{cmd,i}$$
$$\Omega_{cmd,i} = \sqrt{T_{cmd,i}/k_T}$$

Saturação: $\Omega_i \in [0, \Omega_{max}]$, $T_{cmd,i} \in [0, T_{max}]$.

### 6.3 Empuxo e torque estáticos

$$T_i = k_T\,\Omega_i^2$$
$$Q_i = k_Q\,\Omega_i^2 = (0.025\,k_T\,R)\,\Omega_i^2$$

### 6.4 Velocidade induzida (Glauert/Froude axial)

$$v_i = -\frac{V_x}{2} + \sqrt{\left(\frac{V_x}{2}\right)^2 + \frac{T_i}{2\rho A_d}}, \quad V_x = \mathbf V^B \cdot \hat{\mathbf e}_i$$

Em hover ($V_x = 0$): $v_i = \sqrt{T_i/(2\rho A_d)}$.

Para $T_i = mg/4 = 3679$ N: $v_i = 24.3$ m/s.

### 6.5 Velocidade no slipstream (contração de wake)

$$\mathbf V_{slip,s} = \sum_{i \in \mathcal{I}_s} k_s\,v_i\,\hat{\mathbf x}_B$$

onde $\mathcal{I}_s$ é o conjunto de rotores cujo wash banha a superfície $s$. $k_s = 1.0$ (Patterson 2014 conservador).

**Matriz de banho** (4 rotores × 4 superfícies):

|  | wing_R | wing_L | htail | vtail |
|---|---|---|---|---|
| Rotor 1 | ✓ | — | ✓ | — |
| Rotor 2 | ✓ | — | ✓ | ✓ |
| Rotor 3 | — | ✓ | ✓ | — |
| Rotor 4 | — | ✓ | ✓ | ✓ |

### 6.6 Forças e momentos propulsivos

$$\mathbf F^B_{prop} = \sum_{i=1}^{N_r} T_i\,\hat{\mathbf e}_i$$
$$\mathbf M^B_{prop} = \sum_{i=1}^{N_r}\left[ \mathbf r_i \times (T_i\,\hat{\mathbf e}_i) - \mathrm{sign}_i\,Q_i\,\hat{\mathbf e}_i \right]$$

### 6.7 Sanity checks de propulsão

| Verificação | Resultado |
|---|---|
| $\Omega_{max}$: $T = k_T \Omega^2 = T_{max}$ | erro <1% ✓ |
| Trim hover: $\sum T_i \cdot \cos(\text{cant}) = mg$ | erro <5% ✓ |
| Convergência ESC em $5\tau$ | >99% atingido ✓ |
| $v_i$ hover por momentum theory | erro <2% ✓ |
| Reação por torque cancela com 4 rotores pares CCW/CW | $\|M_{react}\|<1$ N·m ✓ |

---

## 7. Arquitetura de Controle

### 7.1 Cascade overview

```
        Trajectory ref (p_ref, v_ref, a_ref, psi_d)
                 ↓
    Outer Loop  PD+FF @ 50 Hz       ω_n = 1.0 rad/s
                 ↓  (f_cmd_NED)
       force_to_attitude  (F → q_d, T_total)
                 ↓
    Inner Loop  SO(3) geometric @ 500 Hz   ω_n = 8.0 rad/s
                 ↓  (M_cmd)
       Control Allocator   pseudoinverse weighted
                 ↓  ([T_1, T_2, T_3, T_4, δe, δa, δr])
       Propulsion + Aerodynamics
                 ↓
       6-DOF Dynamics → RK4 @ 500 Hz
```

### 7.2 Malha externa: PD com feedforward de aceleração

$$\boxed{\mathbf f_{cmd} = -K_p (\mathbf p - \mathbf p_{ref}) - K_v (\mathbf v - \mathbf v_{ref}) + \mathbf a_{ref} - \mathbf g_{NED}}$$

**Tuning** ($\omega_n = 1$ rad/s, $\zeta = 0.85$):
- $K_p = \omega_n^2\,I = 1.0\,I$
- $K_v = 2\zeta\omega_n\,I = 1.7\,I$

**Saturação**: $\|\mathbf f_{cmd}\| \le a_{max} = 2.5 g = 24.5$ m/s² (direção preservada).

**Erro estacionário sob perturbação constante** $d$:
$$\text{err}_{ss} = d/K_p$$
- Hover ($d = 0.16$ m/s²): err ≈ 16 cm
- Cruise ($d = 0.70$ m/s²): err ≈ 70 cm

### 7.3 Mapeamento força → atitude (SE(3))

Dado $\mathbf F_{des,NED}$ e $\psi_d$:

$$\hat{\mathbf x}_B = \frac{\mathbf F_{des}}{\|\mathbf F_{des}\|}$$
$$\hat{\mathbf c}_\psi = [\cos\psi_d; \sin\psi_d; 0]$$
$$\hat{\mathbf y}_B = \frac{\hat{\mathbf c}_\psi \times \hat{\mathbf x}_B}{\|\cdot\|}$$
$$\hat{\mathbf z}_B = \hat{\mathbf x}_B \times \hat{\mathbf y}_B$$
$$R_d = [\hat{\mathbf x}_B,\; \hat{\mathbf y}_B,\; \hat{\mathbf z}_B], \quad \mathbf q_d = \mathrm{from}(R_d)$$

Magnitude do empuxo total comandado: $T_{total} = \|\mathbf F_{des}\| = m \cdot \|\mathbf f_{cmd}\|$.

### 7.4 Malha interna: Controle geométrico SO(3)

**Lei de controle** (Lee-Leok-McClamroch 2010 com integrador):

$$\mathbf M_{cmd} = -k_R\,\mathbf e_R - k_\Omega\,\mathbf e_\Omega - k_I\,\mathbf e_I + \boldsymbol\omega \times J\boldsymbol\omega$$

**Erros**:
$$\mathbf e_R = \frac{1}{2}\mathrm{vee}(R_d^\top R - R^\top R_d)$$
$$\mathbf e_\Omega = \boldsymbol\omega - R^\top R_d \boldsymbol\omega_d$$
$$\dot{\mathbf e}_I = \mathbf e_R, \quad \|\mathbf e_I\| \le I_{max}$$

**Ganhos** (escalados pela inércia, $\omega_n = 8$ rad/s, $\zeta = 0.7$):
- $k_R = J \omega_n^2 = \mathrm{diag}(140\,800,\; 243\,200,\; 352\,000)$
- $k_\Omega = 2\zeta\omega_n J = \mathrm{diag}(24\,640,\; 42\,560,\; 61\,600)$
- $k_I = 0.05 J$, $I_{max} = 15°$

### 7.5 Alocador de controle

**Mapeamento virtual → atuadores** ($6 \times 7$ matriz $B$):

$$\mathbf u = W^{-1} B^\top (B W^{-1} B^\top)^{-1}\,\boldsymbol\nu$$

onde:
$$\boldsymbol\nu = \begin{bmatrix} \|\mathbf F_{cmd}\| \\ 0 \\ 0 \\ \mathbf M_{cmd} \end{bmatrix}, \quad \mathbf u = [T_1, T_2, T_3, T_4, \delta_e, \delta_a, \delta_r]^\top$$

Por rotor: $B_{1:3,i} = \hat{\mathbf e}_i$, $B_{4:6,i} = \mathbf r_i \times \hat{\mathbf e}_i$.

Superfícies: $B_{4:6,5:7}$ escalado por pressão dinâmica $\bar q$ em tempo real.

**Decisão arquitetural crítica**: a primeira coluna de $\boldsymbol\nu$ é a **magnitude** do empuxo ao longo de $+\hat x_B$, **não** $R_{BW}^\top \mathbf F_{cmd,NED}$. Isso é correto porque os rotores fixos só conseguem empurrar em $\hat x_B$ — pedir componentes laterais via alocador (que sairiam via superfícies/diferencial) cria conflito com o SO(3) que está rotacionando o corpo para alinhar $\hat x_B$ com $\mathbf F_{cmd}$.

### 7.6 Limites físicos de autoridade

| Eixo | $J$ (kg·m²) | $M_{max}$ disponível (N·m) | Tempo mínimo para 90° |
|---|---|---|---|
| Roll | 2200 | ~7200 | 7-8 s |
| Pitch | 3800 | ~2700 | 8-10 s |
| Yaw | 5500 | ~90 000 | 1-2 s |

**Implicação**: manobras de transição hover→cruzeiro de 90° em pitch devem ser planejadas com pelo menos 8 segundos para ficar dentro do envelope.

---

## 8. Planejamento de Trajetória

### 8.1 Planicidade diferencial

A aeronave é diferencialmente plana em $\boldsymbol\sigma = [x, y, z, \psi]^\top$. Dada $\boldsymbol\sigma(t) \in C^4$:

$$\mathbf F_{thrust,NED} = m(\ddot{\boldsymbol\sigma}_{xyz} - \mathbf g_{NED})$$
$$\hat{\mathbf x}_B = \mathbf F_{thrust,NED} / \|\cdot\|$$
$$\boldsymbol\omega_d = R_d^\top \frac{\dddot{\boldsymbol\sigma}_{xyz} \times \hat{\mathbf x}_B}{\|\mathbf F\|/m}$$

### 8.2 Polinômios mínimo-snap 8ª ordem

Por segmento $k$ entre waypoints, $C^3$-contínuos:

$$\sigma_d^{(k)}(\tau) = \sum_{j=0}^{7} c_j^{(d,k)}\,\tau^j$$

**Boundary conditions** (8 BCs, 8 coeficientes — sistema bem-determinado):
- $p(0), p'(0), p''(0), p'''(0)$ no início
- $p(T), p'(T), p''(T), p'''(T)$ no fim

### 8.3 Normalização temporal (numérica)

**Problema**: matriz $A$ original tem entradas $T^j$ até $T^7$. Para $T = 116$ s: $T^7 = 3.4 \times 10^{14}$, condicionamento $> 10^{16}$ — coeficientes numericamente inúteis.

**Solução**: resolver em tempo normalizado $u = \tau/T \in [0,1]$. Matriz tem entradas em $\{0, 1, 2, 6, ..., 210\}$, condicionamento $< 100$.

**BCs reescaladas**: $b(2) = T \cdot v_{start}$, $b(3) = T^2 \cdot a_{start}$, $b(4) = T^3 \cdot j_{start}$.

**Avaliação**: $\frac{dp}{d\tau} = \frac{1}{T} \frac{dp}{du}$, e analogamente para derivadas superiores.

### 8.4 Refinamento por algoritmo genético

Custo energético:
$$J_{energy} = \int_0^T \|\mathbf T(t)\|^{1.5}\,dt$$

Decisão: fatores de escala temporal $s_i \in [0.5, 2.0]$ por segmento. GA com população 24, gerações 15.

Redução típica: ~17% (1.252e8 → 1.007e8 nos testes).

### 8.5 Trajetória composável (MissionTrajectory)

Fases concatenáveis:
- **rest_to_rest**: BCs $v=a=j=0$ em ambas pontas
- **rest_to_cruise**: BCs assimétricas (rest início, $v_{cruise}$ fim)
- **cruise_to_rest**: inverso
- **hover**: posição constante, $v=a=0$
- **cruise**: velocidade constante, $a=0$

Continuidade entre fases verificada na construção.

---

## 9. Ambiente e Perturbações

### 9.1 Atmosfera ISA

Modelo padrão para troposfera ($h \le 11$ km):
$$T(h) = T_0 - L h, \quad L = 0.0065 \text{ K/m}$$
$$\rho(h) = \rho_0 \left(\frac{T}{T_0}\right)^{g/(RL) - 1}$$

com $T_0 = 288.15$ K, $\rho_0 = 1.225$ kg/m³, $R = 287.058$ J/(kg·K).

### 9.2 Vento constante + turbulência Dryden (MIL-F-8785C)

**Vento médio**: $\mathbf v_{wind,NED}^{const}$ (default $[3, -1, 0]$ m/s).

**Turbulência Dryden** (filtros formadores em baixa altitude $h < 1000$ ft):

$$H_u(s) = \sigma_u\,\sqrt{\frac{2 L_u}{\pi V}}\,\frac{1}{1 + (L_u/V)s}$$
$$H_{v,w}(s) = \sigma_{v,w}\,\sqrt{\frac{L_{v,w}}{\pi V}}\,\frac{1 + \sqrt{3}(L_{v,w}/V)s}{(1 + (L_{v,w}/V)s)^2}$$

**Comprimentos de escala**:
$$L_w = h_{ft}, \quad L_u = L_v = \frac{h_{ft}}{(0.177 + 0.000823 h_{ft})^{1.2}}$$

**Intensidades** ($W_{20}$ = velocidade do vento a 20 ft):
$$\sigma_w = 0.1 W_{20}, \quad \sigma_u = \sigma_v = \frac{\sigma_w}{(0.177 + 0.000823 h_{ft})^{0.4}}$$

Implementação: discretização Tustin dos filtros, estados iniciais zero, semente RNG fixa para repetibilidade.

---

## 10. Métodos Numéricos

### 10.1 Integração RK4 (passo fixo)

Esquema clássico de 4ª ordem:
$$\mathbf x_{n+1} = \mathbf x_n + \frac{\Delta t}{6}(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4)$$

com:
$$\mathbf k_1 = f(t_n, \mathbf x_n)$$
$$\mathbf k_2 = f(t_n + \Delta t/2, \mathbf x_n + \Delta t/2\,\mathbf k_1)$$
$$\mathbf k_3 = f(t_n + \Delta t/2, \mathbf x_n + \Delta t/2\,\mathbf k_2)$$
$$\mathbf k_4 = f(t_n + \Delta t, \mathbf x_n + \Delta t\,\mathbf k_3)$$

**Renormalização do quatérnio** após cada passo:
$$\mathbf q \leftarrow \mathbf q / \|\mathbf q\|$$

Passo fixo $\Delta t = 2$ ms (500 Hz) — adequado para HIL e compatível com a banda do loop interno.

**Validação empírica**:
- Decaimento $\dot x = -x$: erro a $t=1$ é $3.1 \times 10^{-11}$ (ordem 4 ✓)
- Oscilador harmônico, 1 período: erro de energia $1.1 \times 10^{-14}$ ✓

### 10.2 Solver QP no NMPC (modo opcional)

`fmincon` com algoritmo SQP (`OptimalityTolerance` $10^{-4}$, `MaxIterations` 50) — usado apenas quando `ctrl.outer.type = 'NMPC'` (não-padrão; PD+FF é o default).

### 10.3 Pseudoinversa regularizada (alocador)

$$\mathbf u = W^{-1} B^\top (B W^{-1} B^\top + \epsilon I)^{-1}\,\boldsymbol\nu$$

com $\epsilon = 10^{-3}$ para evitar singularidades quando $\bar q \to 0$ (superfícies sem autoridade).

---

## 11. Tabela Completa de Parâmetros

### 11.1 Constantes físicas

| Símbolo | Valor | Unidade | Descrição |
|---|---|---|---|
| $g$ | 9.80665 | m/s² | Aceleração da gravidade |
| $\rho_0$ | 1.225 | kg/m³ | Densidade ISA SL |
| $T_0$ | 288.15 | K | Temperatura ISA SL |
| $R_{gas}$ | 287.058 | J/(kg·K) | Constante do ar |
| $L$ | 0.0065 | K/m | Gradiente térmico ISA |

### 11.2 Configuração da aeronave

| Símbolo | Valor | Unidade |
|---|---|---|
| $m$ | 1500 | kg |
| $W$ | 14 715 | N |
| $b$ | 11.0 | m |
| $c$ | 1.4 | m |
| $S$ | 15.0 | m² |
| $AR$ | 8.07 | — |
| $e$ | 0.85 | — |
| $C_{L_\alpha}$ | 5.7 | 1/rad |
| $C_{L,max}$ | 1.35 | — |
| $\alpha_{stall}$ | 0.262 (15°) | rad |
| $C_{D_0}$ | 0.022 | — |
| $C_{m_0}$ | 0 | — |
| $C_{m_\alpha}$ | -0.45 | 1/rad |
| $J_{xx}$ | 2200 | kg·m² |
| $J_{yy}$ | 3800 | kg·m² |
| $J_{zz}$ | 5500 | kg·m² |
| $J_{xz}$ | 180 | kg·m² |
| $\mathbf r_{CG}$ | $[0.05; 0; 0]$ | m |

### 11.3 Propulsão (por rotor; 4 rotores idênticos)

| Símbolo | Valor | Unidade |
|---|---|---|
| $R$ | 0.90 | m |
| $A_d$ | 2.545 | m² |
| $N_b$ | 5 | — |
| $c_{blade}$ | 0.10 | m |
| $\theta_{blade}$ | 0.209 (12°) | rad |
| $\sigma$ | 0.177 | — |
| $J_r$ | 0.32 | kg·m² |
| $\Omega_{max}$ | 230 | rad/s |
| $T_{max}$ | 4500 | N |
| $k_T$ | 0.0850 | N·s² |
| $k_Q$ | 0.00191 | N·m·s² |
| $\tau_m$ | 0.060 | s |

### 11.4 Cant dos rotores (varredura paramétrica)

| Rotor | Pitch cant | Yaw cant | $\hat{\mathbf e}_i$ aprox. | Spin |
|---|---|---|---|---|
| 1 | +2° | +1° | [0.999, +0.017, -0.035] | +1 (CCW) |
| 2 | -2° | +1° | [0.999, +0.017, +0.035] | -1 (CW) |
| 3 | +2° | -1° | [0.999, -0.017, -0.035] | -1 (CW) |
| 4 | -2° | -1° | [0.999, -0.017, +0.035] | +1 (CCW) |

### 11.5 Slipstream

| Símbolo | Valor | Unidade |
|---|---|---|
| $k_s$ | 1.0 | — |
| $v_i$ hover | 24.3 | m/s |
| $v_{wake}$ far ($=2v_i$) | 48.6 | m/s |

### 11.6 Superfícies de controle

| Símbolo | Valor |
|---|---|
| $\eta_{elev}$ | 0.55 |
| $\eta_{ail}$ | 0.18 |
| $\eta_{rud}$ | 0.45 |
| $\delta_{e,max}$ | 0.436 (25°) |
| $\delta_{a,max}$ | 0.436 (25°) |
| $\delta_{r,max}$ | 0.436 (25°) |
| $\tau_{servo}$ | 0.040 s |

### 11.7 Controle (ganhos)

**SO(3) Inner** ($\omega_n = 8$ rad/s, $\zeta = 0.7$):
| Ganho | Valor |
|---|---|
| $k_{R}$ | diag(140 800, 243 200, 352 000) |
| $k_{\Omega}$ | diag(24 640, 42 560, 61 600) |
| $k_{I}$ | 0.05 · diag($J$) |
| $I_{max}$ | 0.262 rad (15°) |

**PD+FF Outer** ($\omega_n = 1$ rad/s, $\zeta = 0.85$):
| Ganho | Valor |
|---|---|
| $K_p$ | $I_3$ |
| $K_v$ | $1.7\,I_3$ |
| $a_{max}$ | 24.5 m/s² (2.5g) |

**INDI** (alternativo):
| Ganho | Valor |
|---|---|
| $\omega_{LP,\omega}$ | 50 rad/s |
| $\omega_{LP,u}$ | 50 rad/s |
| $K_{rate}$ | diag(20, 25, 12) |
| $K_{att}$ | diag(10, 12, 6) |

**Alocador**:
| Item | Valor |
|---|---|
| $W$ (effort) | diag(1, 1, 1, 1, 5, 5, 5) |
| Regularização | $\epsilon = 10^{-3}$ |

**Frequências de loop**:
| Loop | Hz |
|---|---|
| Inner SO(3) / RK4 | 500 |
| Outer PDFF | 50 |
| Motor ESC | implicit (1/$\tau_m$ ≈ 17 Hz) |

### 11.8 Ambiente e simulação

| Item | Valor |
|---|---|
| $\Delta t$ simulação | 2 ms |
| $t_{final}$ (missão) | 600 s |
| Vento médio NED | [3, -1, 0] m/s |
| $W_{20}$ Dryden | 7.7 m/s |
| Dryden enable | true |
| RNG seed | 42 |

---

## 12. Validação e Limites

### 12.1 Sanity checks físicos (`sanity_check_model.m`)

| # | Métrica | Valor | Status |
|---|---|---|---|
| 1 | T/W ratio | 1.22 | ✅ OK |
| 2 | Hover throttle | 81.7% | ✅ OK |
| 3 | Wing loading | 100 kg/m² | ✅ OK |
| 4 | Disk loading | 147 kg/m² | ✅ OK |
| 5 | Hover power | 510 kW | ✅ OK |
| 6 | Power loading | 2.94 kg/kW | ✅ OK |
| 7 | $V_{stall}$ | 34.4 m/s | ✅ OK |
| 8 | $V_{LDmax}$ | 48.2 m/s | ✅ OK |
| 9 | $\alpha_{cruise}$ | 8.0° | ✅ OK |
| 10 | $L/D_{cruise}$ | 15.5 | ✅ OK |
| 11 | Bandwidth inner/motor | 3.0× | ✅ OK |
| 12 | Bandwidth outer/inner | 8.0× | ✅ OK |
| 13 | $J_{yy}/J_{xx}$ | 1.73 | ✅ OK |
| 14 | $J_{xz}/J_{yy}$ | 0.047 | ✅ OK |
| 15 | Pitch authority | 2700 N·m | ⚠️ Limite físico (não é bug) |
| 16 | Crossover aero/rotor | 40 m/s | ✅ OK |
| 17 | Queda livre $a = g$ | $10^{-12}$ | ✅ OK |
| 18 | Trim hover $T = mg$ | $10^{-9}$ | ✅ OK |

**Pitch authority**: 2700 N·m disponíveis vs ~18 700 N·m necessários para pitch 90°/1s. **Implica tempo mínimo de manobra de transição = 8 s** (característica de design do veículo, não bug).

### 12.2 Suíte de testes automatizados

| Etapa | Suíte | Testes | Status |
|---|---|---|---|
| 1 | quat_utils, so3_utils, dinâmica | 27 | ✅ PASS |
| 2 | aerodinâmica, propulsão | 17 | ✅ PASS |
| 3 | SO(3), INDI, alocador, force_to_attitude | 33 | ✅ PASS |
| 4 | Flatness, trajectory optimizer | 15 | ✅ PASS |
| 5 | RK4, Dryden, atmosfera | 11 | ✅ PASS |
| 6 | Integração hover (10 s) | 6 | ✅ PASS (drift 0.18 m) |
| 7 | Integração missão dinâmica | 7 | ⏳ Em validação com PD+FF |
| **Total** | **17 suítes** | **118** | **111 PASS** |

### 12.3 Faixa de validade

| Aspecto | Limite |
|---|---|
| Mach | ≤ 0.3 (incompressível) |
| Altitude | ≤ 11 km (troposfera) |
| Ângulo de ataque | ±90° (Viterna-Corrigan) |
| Velocidade vertical | ≤ 30 m/s (validade do momentum theory) |
| Transição pitch | ≥ 8 s (limite de autoridade) |
| Vento turbulento | $W_{20}$ ≤ 15 m/s (categoria leve/moderado) |

---

## 13. Referências

[1] Stevens, B.L., Lewis, F.L., Johnson, E.N. (2015). *Aircraft Control and Simulation*, 3rd ed., Wiley.

[2] Etkin, B., Reid, L.D. (1996). *Dynamics of Flight: Stability and Control*, 3rd ed., Wiley.

[3] Lustosa, L.R., Defaÿ, F., Moschetta, J.-M. (2017). *Global Singularity-Free Aerodynamic Model for Algorithmic Flight Control of Tail Sitters*. J. Guidance Control Dynamics 40(7).

[4] Markley, F.L., Crassidis, J.L. (2014). *Fundamentals of Spacecraft Attitude Determination and Control*. Springer.

[5] Viterna, L.A., Corrigan, R.D. (1981). *Fixed Pitch Rotor Performance of Large Horizontal Axis Wind Turbines*. NASA CP-2230.

[6] Anderson, J.D. (1991). *Fundamentals of Aerodynamics*, 2nd ed., McGraw-Hill.

[7] Patterson, M.D., Borer, N.K., German, B. (2014). *Wing-Propeller Interaction: Validation of Computational Methods*. NASA TM.

[8] Borer, N.K., Patterson, M.D., et al. (2016). *Design and Performance of the NASA SCEPTOR Distributed Electric Propulsion Flight Demonstrator*. AIAA Aviation Forum.

[9] Leishman, J.G. (2006). *Principles of Helicopter Aerodynamics*, 2nd ed., Cambridge.

[10] Drela, M. (2014). *Flight Vehicle Aerodynamics*. MIT Press.

[11] Bristeau, P.-J., Martin, P., Salaün, E., Petit, N. (2009). *The role of propeller aerodynamics in the model of a quadrotor UAV*. ECC.

[12] Lee, T., Leok, M., McClamroch, N.H. (2010). *Geometric tracking control of a quadrotor UAV on SE(3)*. IEEE CDC.

[13] Mellinger, D., Kumar, V. (2011). *Minimum snap trajectory generation and control for quadrotors*. IEEE ICRA.

[14] Smeur, E.J.J., Chu, Q., de Croon, G.C.H.E. (2016). *Adaptive Incremental Nonlinear Dynamic Inversion for Attitude Control of MAVs*. J. Guidance Control Dynamics 39(3).

[15] Mellinger, D., Mueller, M.W., Kumar, V. (2012). *Trajectory generation and control for precise aggressive maneuvers with quadrotors*. Int. J. Robotics Research 31(5).

[16] Bouabdallah, S., Murrieri, P., Siegwart, R. (2004). *Design and control of an indoor micro quadrotor*. IEEE ICRA.

[17] U.S. Department of Defense (1980). *Military Specification: Flying Qualities of Piloted Airplanes*, MIL-F-8785C.

[18] U.S. National Oceanic and Atmospheric Administration (1976). *U.S. Standard Atmosphere*. NOAA-S/T 76-1562.

[19] Joby Aviation. *Joby S4 Aircraft Specifications*. https://www.jobyaviation.com (acesso: 2026).

[20] Beta Technologies. *Beta Alia 250 Specifications*. https://www.beta.team (acesso: 2026).

[21] Verling, S., Weibel, B., Boosfeld, M., Alexis, K., Burri, M., Siegwart, R. (2016). *Full attitude control of a VTOL tailsitter UAV*. IEEE ICRA.

---

## Apêndice A: Estrutura de Software (~4500 linhas)

```
evtol_simulation/
├── main.m (75 linhas)
├── config/
│   ├── aircraft_config.m       (~140 linhas)
│   ├── controller_config.m     (~70 linhas)
│   └── simulation_config.m     (~50 linhas)
├── core/
│   ├── Aircraft.m              (~130 linhas)
│   ├── quat_utils.m            (~150 linhas)
│   └── so3_utils.m             (~70 linhas)
├── aerodynamics/
│   ├── Aerodynamics.m          (~200 linhas)
│   └── Propulsion.m            (~140 linhas)
├── control/
│   ├── AttitudeControllerSO3.m (~80 linhas)
│   ├── AttitudeControllerINDI.m (~110 linhas)
│   ├── PositionControllerPDFF.m (~70 linhas)   [DEFAULT]
│   ├── PositionControllerNMPC.m (~180 linhas)
│   ├── ControlAllocator.m      (~110 linhas)
│   └── force_to_attitude.m     (~100 linhas)
├── trajectory/
│   ├── DifferentialFlatness.m  (~130 linhas)
│   ├── MissionTrajectory.m     (~200 linhas)
│   └── TrajectoryOptimizer.m   (~120 linhas)
├── environment/
│   ├── DrydenWind.m            (~110 linhas)
│   └── atmosphere_isa.m        (~25 linhas)
├── simulation/
│   ├── run_simulation.m        (~210 linhas)
│   ├── rk4_step.m              (~35 linhas)
│   ├── build_mission.m         (~50 linhas)
│   └── build_realistic_mission.m (~110 linhas)
├── visualization/
│   └── plot_telemetry.m        (~110 linhas)
├── tests/                      (~1300 linhas)
└── docs/                       (~2000 linhas Markdown)
```

---

*Fim do Relatório Técnico v1.0*
