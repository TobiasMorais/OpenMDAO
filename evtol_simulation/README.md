# Plataforma de Simulação Virtual para Tailsitter eVTOL Biplace 1500 kg

Framework MATLAB orientado a objetos para simulação 6-DOF de aeronaves eVTOL não-convencionais (foco em **tailsitter biplace 1500 kg, 4 rotores fixos com cant configurável**), cobrindo todo o envelope de voo com ênfase na transição pós-estol hover→cruzeiro.

Implementa as cinco etapas modulares definidas no relatório técnico:

| Etapa | Conteúdo |
|---|---|
| 1 | 6-DOF Newton-Euler com quatérnios + tensor com $J_{xz}$ |
| 2 | Aerodinâmica Viterna-Corrigan 360° + BEMT + slipstream + giroscópio rotor |
| 3 | Cascata: Geométrico SO(3) **OU** INDI (interno) + NMPC (externo) + alocador |
| 4 | Trajetória ótima: Planicidade Diferencial (min-snap) + GA refinador energético |
| 5 | Loop RK4 fixed-step 500 Hz + turbulência Dryden + plots Grau Aeroespacial |

Detalhes matemáticos completos: [`docs/EQUATIONS.md`](docs/EQUATIONS.md)

## Estrutura

```
evtol_simulation/
├── main.m                          # entry point (run this)
├── config/
│   ├── aircraft_config.m           # 1500 kg tailsitter parâmetros físicos
│   ├── controller_config.m         # ganhos SO(3)/INDI/NMPC, alocador
│   └── simulation_config.m         # RK4, Dryden, cenário
├── core/
│   ├── Aircraft.m                  # classe principal (6-DOF dynamics)
│   ├── quat_utils.m                # quatérnios Hamilton (mul, conj, kinematic, errMul, exp/log)
│   └── so3_utils.m                 # álgebra de Lie so(3) (hat, vee, errMat, exp, log)
├── aerodynamics/
│   ├── Aerodynamics.m              # strip theory + Viterna-Corrigan 360°
│   └── Propulsion.m                # 4 rotores: 1ª ordem ESC, BEMT, slipstream, M_gyro
├── control/
│   ├── AttitudeControllerSO3.m     # geométrico Lee-Leok-McClamroch + integrador so(3)
│   ├── AttitudeControllerINDI.m    # incremental dynamic inversion (sensor-driven)
│   ├── PositionControllerNMPC.m    # SQP fmincon, point-mass NED, tilt envelope
│   ├── ControlAllocator.m          # pseudoinversa ponderada 6→7 atuadores
│   └── force_to_attitude.m         # mapeamento F_NED → (q_d, T_total, omega_d)
├── trajectory/
│   ├── DifferentialFlatness.m      # min-snap polinômios em (x,y,z,psi)
│   └── TrajectoryOptimizer.m       # baseline flat + GA energético (∫T^1.5 dt)
├── environment/
│   ├── DrydenWind.m                # MIL-F-8785C low-altitude
│   └── atmosphere_isa.m            # densidade ISA troposférica
├── simulation/
│   ├── rk4_step.m                  # RK4 com renormalização de quatérnio
│   ├── build_mission.m             # waypoints (full_mission/transition_only/hover_disturbed)
│   └── run_simulation.m            # closed-loop cascata
├── visualization/
│   └── plot_telemetry.m            # 3D, NED, PQR, motor diff, pitch, erro
└── docs/
    └── EQUATIONS.md                # equacionamento completo em LaTeX
```

## Uso rápido

### Passo 1 — Rodar a suíte de testes ANTES de qualquer simulação

```matlab
cd evtol_simulation
addpath('tests')
tests/run_all_tests
```

A suíte (15 testes, 92 verificações) valida cada componente isoladamente. **Toda integração deve esperar 15/15 PASS**. Detalhes: [`docs/TEST_PROCEDURE.md`](docs/TEST_PROCEDURE.md).

### Passo 2 — Rodar a missão completa

```matlab
main
```

Editar `main.m` para alternar:

```matlab
USE_GA_REFINER = true;       % otimização energética via GA (mais lento)
INNER_LAW      = 'INDI';     % 'SO3' ou 'INDI'
sim.scenario   = 'full_mission';   % 'transition_only' | 'hover_disturbed'
```

## Cenários default

- **`full_mission`**: takeoff vertical → climb → transição → cruzeiro 1.2 km → desaceleração → pouso vertical (~120 s)
- **`transition_only`**: foco na manobra crítica hover→cruise (~15 s)
- **`hover_disturbed`**: pairado estacionário com Dryden ativo, valida robustez

## Cant dos rotores (parâmetro de varredura)

Os 4 rotores possuem ângulos $(\alpha_p^i, \alpha_y^i)$ configuráveis em `aircraft_config.m`:

```matlab
ac.rotor.pitch_cant = deg2rad([+2.0; -2.0; +2.0; -2.0]);
ac.rotor.yaw_cant   = deg2rad([+1.0; +1.0; -1.0; -1.0]);
```

Variar entre execuções para estudar trade-offs de design (autoridade vs perda de eficiência).

## Validação Grau Aeroespacial (telemetrias)

Conforme requisito original e PDF [seção HIL]:
- Posição NED, Taxas PQR, Esforço Diferencial de Motores, Ângulo de Arfagem, Erro de Rastreamento — todos plotados.
- Energy proxy $\int T^{1.5}\,dt$ impresso ao final.
- Norma do quaternion garantida por renormalização pós-RK4.

## Dependências MATLAB

- **Base**: MATLAB R2020a+
- **Optimization Toolbox** (fmincon para NMPC) — *recomendado*
- **Global Optimization Toolbox** (ga para refinador energético) — *opcional*

Sem fmincon, NMPC cai automaticamente em fallback PD interno.
Sem ga, o framework usa apenas trajetória baseline de planicidade.

## Próximos passos (roadmap)

- VLM de baixa ordem para acoplamento prop-asa com swirl
- Sensores IMU/GPS simulados com viés e ruído (já parcialmente em config)
- Interface UDP/TCP para Hardware-in-the-Loop
- Modos de falha de rotor (perda de empuxo) para validação de robustez
- Conversão para Simulink/S-Function para deploy em RTOS

## Referências (do PDF)

Lee-Leok-McClamroch SO(3) [27], Smeur INDI [67-68], Mellinger flatness [73],
Viterna-Corrigan post-stall [56], Patterson slipstream [53], MIL-F-8785C Dryden,
Joby/Beta scaling refs [2,17,21,46,55].
