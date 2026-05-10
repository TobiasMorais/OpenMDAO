# Procedimento Completo de Simulação — Validação Sistemática

**Objetivo**: validar tudo (matemática, física, código, integração) **em camadas**, do menor ao maior. Se uma camada falhar, parar e diagnosticar antes de prosseguir.

```
        Camada 5: Validação física pós-simulação (validate_physics)
           ↑
        Camada 4: Missão completa visual (main)
           ↑
        Camada 3: Integração (test_integration_hover, test_integration_realistic)
           ↑
        Camada 2: Sanidade física do modelo (sanity_check_model)
           ↑
        Camada 1: Testes unitários por módulo (run_all_tests)
```

## Pré-requisitos

```matlab
cd evtol_simulation
addpath(genpath(pwd));
```

## CAMADA 1 — Testes unitários

**Comando**:
```matlab
results = run_all_tests;
```

**Critério de aprovação**: `results.ready_for_integration == 1`, ou seja **16/16 PASS, 111/111 testes**.

Suítes:
| Etapa | Suíte | Verifica |
|---|---|---|
| 1 | `test_quat_utils` | Algebra de quatérnios (Hamilton mul, conj, exp/log, kinematic) |
| 1 | `test_so3_utils` | hat/vee, errMat, Rodrigues, det |
| 1 | `test_aircraft_dynamics` | Newton-Euler, queda livre, hover trim, J_xz |
| 2 | `test_aerodynamics` | Viterna 360°, regime linear, slipstream |
| 2 | `test_propulsion` | T = kT·Ω², ESC τ, momentum theory v_i |
| 3 | `test_so3_controller` | Bandwidth, damping, settling time |
| 3 | `test_indi_controller` | Filtros, condicionamento G |
| 3 | `test_nmpc` | Hover, step, restrição de tilt |
| 3 | `test_allocator` | Pseudoinversa, momentos, saturação |
| 3 | `test_force_to_attitude` | F→qd, det=+1, continuidade, heading |
| 4 | `test_flatness` | BCs polinômio, derivadas, continuidade |
| 4 | `test_trajectory_opt` | Baseline, GA, factibilidade |
| 5 | `test_rk4` | Ordem 4, q-renorm, energia |
| 5 | `test_dryden` | Variância, reset, zero-mean |
| 5 | `test_atmosphere` | ISA, monotonicidade |
| 6 | `test_integration_hover` | Cascata em hover (corpo não flipa, drift < 5 m) |
| 7 | `test_integration_realistic` | Missão completa por fases |

**Se falhar**: relatório imprime exatamente quais asserções falharam, com tolerâncias e valores. Cada teste tem comentário inline explicando a física.

## CAMADA 2 — Sanidade física do modelo

**Comando**:
```matlab
results = sanity_check_model;
```

**Critério**: 22+ marcadores `[OK]`, máximo 1-2 `[WARN]`, **0 `[FAIL]`**.

Verifica 23 indicadores físicos derivados de `aircraft_config.m`:

| Categoria | Métricas |
|---|---|
| Sizing | T/W, hover throttle, wing loading, disk loading |
| Performance | V_stall, V_cruise, CL_cruise, α_cruise, L/D, P_hover, P_cruise |
| Dinâmica | Separação de bandwidth (motor/inner/outer) |
| Inércia | $J_{yy}/J_{xx}$, $J_{xz}/J_{yy}$ |
| Autoridade | Pitch moment max vs requerido, crossover aero/rotor |
| Conservação | Queda livre $a=g$, hover trim $T=mg$ |

**Se algum FAIL**: parâmetro físico inconsistente (ex: aeronave que não consegue pairar, T_max < m·g/4). Verificar `aircraft_config.m`.

## CAMADA 3 — Integração da cascata

Já incluída no `run_all_tests` (Stages 6 e 7), mas pode ser rodada isoladamente:

```matlab
test_integration_hover         % 10 s hover, drift < 5 m, pitch < 15° dev
test_integration_realistic     % missão completa, ~250 s
```

**Comportamento esperado**:
- Hover: drift < 1 m, pitch dev < 5°, motores não saturam
- Missão realista:
  - Fase 1 (climb 20m): erro < 20 m
  - Fase 2 (hover 20s): erro < 5 m
  - Fase 3 (climb-trans 1000m): erro < 200 m durante transiente
  - Fase 4 (cruise): erro < 100 m após assentamento
  - Pitch em cruise: ~ α_LDmax (6.9°)

## CAMADA 4 — Missão completa visual

**Comando**:
```matlab
main
```

**Saída esperada**:
1. **Console**: resumo da missão com V_LDmax, durações por fase, total
2. **Console**: relatório de tracking ao final (mean, max, final error)
3. **Console**: validação física automática `validate_physics`
4. **Plots** (8 figuras):
   - 3D Trajectory: caminho próximo do planejado
   - Position NED: monotônico até cruise, depois constant ALT
   - Velocity NED: V_N → 48 m/s em cruise
   - Body rates p,q,r: picos em transição, calmo em cruise
   - Pitch θ: 90° → 7° em cruise, → 90° de volta em hover (se phases extras)
   - Per-rotor thrust: variação suave, sem saturação prolongada
   - Differential motor effort: pequeno em cruise estável
   - Tracking error: bounded ao longo da missão

**Indicadores de problema**:
- Pitch oscilando ±90°: cascata instável
- Rotores em 4500 N o tempo todo: trajetória demanda mais que envelope
- Tracking error crescendo monotonicamente: divergência

## CAMADA 5 — Validação física pós-simulação

Roda automaticamente dentro de `main`, mas pode ser feita standalone:

```matlab
validate_physics(log, ac, mission_info)
```

Compara saídas da simulação com **predições analíticas**:

| ID | Verificação | Predição |
|---|---|---|
| V1 | Hover thrust trim | $T_{each} \approx mg/4 = 3679$ N (erro < 5%) |
| V2 | Quaternion norm | $||q|| - 1 < 10^{-3}$ sempre |
| V3 | Cruise drag balance | $D_{pred} \approx T_{horiz}$ a partir do polar |
| V4 | Energy budget | $W_{thrust} \ge \Delta PE \cdot 0.5$ |
| V5 | Cruise pitch | $\theta \approx \alpha_{LDmax} = 6.9°$ |
| V6 | Tracking statistics | mean, max, final error |
| V7 | Motor utilization | avg < 80%, sem saturação prolongada |
| V8 | Final attitude | euler angles físicamente plausíveis |

## Estrutura de Arquivos Reference

```
evtol_simulation/
├── main.m                         # entry point (Camadas 4-5)
├── docs/
│   ├── EQUATIONS.md               # equações de cada etapa
│   ├── MISSION_DESIGN.md          # detalhe missão realista
│   ├── SANITY_CHECKS.md           # documentação Camada 2
│   ├── TEST_PROCEDURE.md          # documentação Camada 1
│   └── RUN_PROCEDURE.md           # este arquivo
├── tests/
│   ├── run_all_tests.m            # Camada 1
│   ├── sanity_check_model.m       # Camada 2
│   ├── validate_physics.m         # Camada 5
│   ├── test_integration_hover.m   # Camada 3a
│   ├── test_integration_realistic.m # Camada 3b
│   └── ... (15 outros testes unitários)
├── core/                          # Aircraft, quat_utils, so3_utils
├── aerodynamics/                  # Aerodynamics, Propulsion
├── control/                       # SO3, INDI, NMPC, Allocator, force_to_attitude
├── trajectory/                    # DifferentialFlatness, MissionTrajectory, TrajectoryOptimizer
├── environment/                   # DrydenWind, atmosphere_isa
├── simulation/                    # run_simulation, rk4_step, build_*
└── visualization/                 # plot_telemetry
```

## Sequência Recomendada de Execução

**Primeira vez** (validação completa):
```matlab
% 1. Sanidade física
sanity_check_model;
% 2. Testes unitários (espera 16/16)
results = run_all_tests;
% 3. Missão completa
main;
% 4. Inspecionar plots
% 5. Resultados validate_physics
```

**Iteração** (após mudar parâmetros):
```matlab
sanity_check_model;     % 5 segundos
main;                    % ~3 minutos
```

**Depuração** (após edição em código):
```matlab
test_quat_utils;        % ou outro teste específico do módulo modificado
test_integration_hover;  % se mudou cascata
```

## Critérios de Aceitação Final ("Pronto para entrega")

✅ Camada 1: 16/16 PASS, 111/111 tests
✅ Camada 2: 22+ OK, 0 FAIL
✅ Camada 3: hover error < 1 m, pitch dev < 5°
✅ Camada 4: missão completa visualmente correta, motores não saturados
✅ Camada 5: V1-V8 todos OK ou INFO

Quando todas as camadas passarem, o sistema está **pronto para uso operacional** — comparar com aeronaves de referência (Joby S4, Beta Alia), variar parâmetros (cant angles, mass), ou estender com novas missões.
