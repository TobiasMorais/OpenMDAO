# Equacionamento Matemático Completo
## Plataforma de Simulação Tailsitter eVTOL Biplace 1500 kg

Todos os símbolos e convenções seguem o relatório técnico de referência (PDF) e os papers citados nas referências [1]–[81].

---

## ETAPA 1 — Cinemática e Dinâmica de Corpo Rígido (6-DOF)

### 1.1 Sistemas de referência

- **Inercial NED** (Frame I): $X_I$=Norte, $Y_I$=Leste, $Z_I$=Down (paralelo a $\vec g$).
- **Corpo** (Frame B): $X_B$=longitudinal (nariz), $Y_B$=ponta direita, $Z_B$=ventral.

### 1.2 Quatérnios unitários (Hamilton, escalar-primeiro)

Estado de atitude $q=[q_0,q_1,q_2,q_3]^\top$ com $\|q\|=1$. Matriz de rotação corpo→mundo:

$$
R(q) = (q_0^2 - \mathbf{q}_v^\top \mathbf{q}_v)\,I + 2\,\mathbf{q}_v\mathbf{q}_v^\top - 2\,q_0\,[\mathbf{q}_v]_\times
$$

Cinemática quaternion (sem singularidades):

$$
\dot q = \tfrac{1}{2}\,\Omega(\boldsymbol\omega^B)\,q,\quad
\Omega(\omega) = \begin{bmatrix}
0 & -\omega_x & -\omega_y & -\omega_z\\
\omega_x & 0 & \omega_z & -\omega_y\\
\omega_y & -\omega_z & 0 & \omega_x\\
\omega_z & \omega_y & -\omega_x & 0
\end{bmatrix}
$$

**Erro multiplicativo** (anti-unwinding):
$$q_e = q_d^{-1} \otimes q,\qquad \text{se } q_{e,0}<0\;\Rightarrow\;q_e\leftarrow -q_e$$

### 1.3 Equações de Newton-Euler (6-DOF)

**Translação no NED** (forma usada no integrador):
$$
\dot{\mathbf p}_{NED} = \mathbf v_{NED},\qquad
\dot{\mathbf v}_{NED} = \frac{1}{m}\,R(q)\,\big(\mathbf F^B_{aero} + \mathbf F^B_{prop}\big) + \mathbf g_{NED}
$$

com $\mathbf g_{NED}=[0,0,g]^\top$, $g=9.80665$ m/s².

**Rotação no corpo** (com $J_{xz}$ acoplado):
$$
J\,\dot{\boldsymbol\omega}^B + \boldsymbol\omega^B\times(J\boldsymbol\omega^B) + \mathbf M^B_{gyro}
=\mathbf M^B_{aero}+\mathbf M^B_{prop}
$$

Tensor de inércia (assimetria fuselagem/cauda):
$$
J=\begin{bmatrix} J_{xx} & 0 & -J_{xz}\\ 0 & J_{yy} & 0\\ -J_{xz} & 0 & J_{zz}\end{bmatrix}
$$

### 1.4 Torque giroscópico de rotores

$$
\mathbf M^B_{gyro} = \boldsymbol\omega^B\times \sum_{i=1}^{N_r} J_{r}\,\Omega_i^{(s)}\,\hat{\mathbf e}_i
$$
onde $\Omega_i^{(s)} = \mathrm{sign}_i\,\Omega_i$ e $\hat{\mathbf e}_i$ é o eixo de empuxo do rotor $i$ no corpo (com cant aplicado).

---

## ETAPA 2 — Modelo Aero-Propulsivo

### 2.1 Aerodinâmica 360° (Viterna-Corrigan)

Regime linear ($|\alpha|\le\alpha_s$):
$$
C_L = C_{L_\alpha}(\alpha-\alpha_0),\qquad
C_D = C_{D_0} + \frac{C_L^2}{\pi A R\,e}
$$

Pós-estol até $\pm 90°$:
$$
C_L(\alpha)=A_1\sin(2\alpha)+A_2\,\frac{\cos^2\alpha}{\sin\alpha},\quad
C_D(\alpha)=B_1\sin^2\alpha+B_2\cos\alpha
$$
com $C_{D,\max}\approx 1.11+0.018\,AR$ e
$$
A_1=\tfrac{C_{D,\max}}{2},\;
A_2=(C_{L,s}-C_{D,\max}\sin\alpha_s\cos\alpha_s)\,\tfrac{\sin\alpha_s}{\cos^2\alpha_s},\;
B_1=C_{D,\max},\;
B_2=\tfrac{C_{D,s}-C_{D,\max}\sin^2\alpha_s}{\cos\alpha_s}
$$

### 2.2 Propulsão (atuador 1ª ordem + BEMT)

Dinâmica do ESC/motor:
$$\tau_m\dot\Omega_i + \Omega_i = \Omega_{cmd,i}$$

Estática:
$$T_i = k_T\,\Omega_i^2,\qquad Q_i = k_Q\,\Omega_i^2$$

Velocidade induzida (Glauert/Froude axial):
$$
v_i = -\frac{V_x}{2} + \sqrt{\frac{V_x^2}{4} + \frac{T_i}{2\rho A_d}},\quad V_x=\mathbf V^B\!\cdot\hat{\mathbf e}_i
$$

Velocidade no slipstream contraído:
$$\mathbf V_{\text{slip},s}=\sum_{i\in\mathcal{I}_s} k_s\,v_i\,\hat{\mathbf e}_i$$

### 2.3 Forças e momentos por strip

Para cada strip $k$ (asa-R, asa-L, htail, vtail):
$$
\mathbf V_{loc,k}=\mathbf V^B+\boldsymbol\omega^B\times\mathbf r_k+\mathbf V_{\text{slip},k}
$$
$$
\bar q_k=\tfrac{1}{2}\rho\|\mathbf V_{loc,k}\|^2,\quad
L_k=\bar q_k S_k C_L(\alpha_k),\quad D_k=\bar q_k S_k C_D(\alpha_k)
$$
Soma:
$$\mathbf F^B_{aero}=\sum_k \mathbf F_k,\qquad \mathbf M^B_{aero}=\sum_k \mathbf r_k\times\mathbf F_k$$

### 2.4 Eixo de empuxo dos rotores com cant

Para cant $(\alpha_p^i,\alpha_y^i)$:
$$\hat{\mathbf e}_i = R_z(\alpha_y^i)\,R_y(\alpha_p^i)\,\hat{\mathbf x}_B$$

---

## ETAPA 3 — Arquitetura de Controle (GNC)

### 3.1 Malha externa: Nonlinear MPC

Modelo preditor (point-mass NED):
$$
\mathbf p_{k+1}=\mathbf p_k+\Delta t\,\mathbf v_k,\qquad
\mathbf v_{k+1}=\mathbf v_k+\Delta t\,(\mathbf f_{cmd,k}+\mathbf g_{NED})
$$

Custo:
$$
J=\sum_{k=0}^{N-1}\Big[(\mathbf p_k-\mathbf p_{r,k})^\top Q_p(\cdot)+(\mathbf v_k-\mathbf v_{r,k})^\top Q_v(\cdot)+\mathbf f_k^\top R\,\mathbf f_k\Big]+\mathbf e_N^\top Q_f \mathbf e_N
$$
sujeito a $|\mathbf f_k|\le a_{\max}$ e $\arctan(\|f_{xy}\|/f_z)\le \mathrm{tilt}_{\max}$.

### 3.2 Malha interna A: Geométrica em SO(3)

Erros (Lee, Leok, McClamroch 2010):
$$
\mathbf e_R = \tfrac{1}{2}\,\mathrm{vee}(R_d^\top R - R^\top R_d),\quad
\mathbf e_\Omega = \boldsymbol\omega - R^\top R_d \boldsymbol\omega_d
$$

Lei de controle com integrador em $\mathfrak{so}(3)$:
$$
\boxed{\,\mathbf M_{cmd}= -k_R\,\mathbf e_R - k_\Omega\,\mathbf e_\Omega - k_I\,\mathbf e_{I} + \boldsymbol\omega\times J\boldsymbol\omega - J(\widehat{\boldsymbol\omega}\,R^\top R_d \boldsymbol\omega_d - R^\top R_d\dot{\boldsymbol\omega}_d)\,}
$$

### 3.3 Malha interna B: INDI (Smeur 2016)

Dinâmica:
$$\dot{\boldsymbol\omega}=f(x)+G\,\mathbf u,\quad G\in\mathbb R^{3\times m}$$

Comando incremental:
$$
\Delta \mathbf u = G^{+}\big(\boldsymbol\nu - \dot{\boldsymbol\omega}_{\text{filt}}\big),\qquad
\mathbf u_{cmd} = \mathbf u_{\text{filt}} + \Delta\mathbf u
$$
com lei sintética $\boldsymbol\nu = K_\Omega(\boldsymbol\omega_{des}-\boldsymbol\omega)+\dot{\boldsymbol\omega}_{des}$ e $\boldsymbol\omega_{des}=-K_\phi\,\mathrm{vec}(q_e)$.

### 3.4 Alocação de controle

$$\boldsymbol\nu = [\mathbf F^B; \mathbf M^B] = B(q,\rho,V)\,\mathbf u,\qquad \mathbf u=[T_1,...,T_4,\delta_e,\delta_a,\delta_r]^\top$$

$$\mathbf u^* = W^{-1}B^\top(BW^{-1}B^\top)^{-1}\boldsymbol\nu \;\;\text{(pseudoinversa ponderada + saturação)}$$

Por rotor: $B_{1:3,i}=\hat{\mathbf e}_i$, $B_{4:6,i}=\mathbf r_i\times\hat{\mathbf e}_i$.

---

## ETAPA 4 — Otimização de Trajetória

### 4.1 Planicidade Diferencial

Saídas planas: $\boldsymbol\sigma=[x,y,z,\psi]^\top$. A partir de $\sigma(t)$ até a 4ª derivada:
$$
\mathbf F^{NED}_{thrust} = m(\ddot{\boldsymbol\sigma}_{xyz}-\mathbf g_{NED}),\qquad
\hat{\mathbf x}_B = \frac{\mathbf F^{NED}_{thrust}}{\|\cdot\|}
$$
$$
\boldsymbol\omega_d = R_d^\top\,\frac{\dddot{\boldsymbol\sigma}_{xyz}\times\hat{\mathbf x}_B}{\|\mathbf F\|/m}
$$

Polinômios mínimos-snap por segmento (8ª ordem, $C^3$ contínuos):
$$
\sigma_d^{(k)}(\tau)=\sum_{j=0}^{7} c_j^{(d,k)}\tau^j,\quad d\in\{x,y,z,\psi\}
$$

### 4.2 Custo energético + GA

$$
\boxed{\,J_{energy}=\int_0^{T} \|\mathbf T(t)\|^{1.5}\,dt\,}
$$

GA otimiza alocação temporal $\mathbf t_{seg}\to \mathbf s\odot\mathbf t_{seg,0}$ com $s_i\in[0.5,2]$, sob restrições $T_i\le T_{\max}$, tilt $\le 89°$.

---

## ETAPA 5 — Ambiente de Simulação

### 5.1 Integração RK4 (passo fixo)

$$
x_{n+1}=x_n+\tfrac{\Delta t}{6}(k_1+2k_2+2k_3+k_4)
$$
com renormalização do quaternion após cada passo: $q\leftarrow q/\|q\|$.

### 5.2 Turbulência Dryden (MIL-F-8785C)

Filtros formadores (baixa altitude $h<1000$ ft):
$$
H_u(s)=\sigma_u\sqrt{\frac{2L_u}{\pi V}}\,\frac{1}{1+\frac{L_u}{V}s},\quad
H_{v,w}(s)=\sigma_{v,w}\sqrt{\frac{L_{v,w}}{\pi V}}\,\frac{1+\sqrt 3\,\frac{L_{v,w}}{V}s}{(1+\frac{L_{v,w}}{V}s)^2}
$$
$$
L_w=h_{ft},\quad L_u=L_v=\frac{h_{ft}}{(0.177+0.000823 h_{ft})^{1.2}}
$$
$$
\sigma_w=0.1\,W_{20},\quad \sigma_u=\sigma_v=\frac{\sigma_w}{(0.177+0.000823 h_{ft})^{0.4}}
$$

### 5.3 Telemetrias logadas (Grau Aeroespacial)

| Sinal | Símbolo | Uso |
|---|---|---|
| Posição NED | $\mathbf p_{NED}$ | Trajetória 3D |
| Taxas PQR | $\boldsymbol\omega^B$ | Estabilidade rotacional |
| Esforço motor | $T_i$ | Saturação, falha |
| Diferencial motor | $T_i-\bar T$ | Acoplamento giroscópico |
| Pitch | $\theta$ | Validação 90°→0° |
| Erro tracking | $\|\mathbf p-\mathbf p_r\|$ | Performance MPC |

### 5.4 Critérios de validação HIL/DO-178C

- $\sup\|q\|-1| < 10^{-6}$ (norma do quaternion)
- LCO ausente em hover com Dryden ativo (10 min sim)
- Erro de tracking < 2 m em cruzeiro estabilizado
- Margem MUAD respeitada na transição (ângulo equivalente)
