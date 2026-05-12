# EVTOL-Tailsitter

**Simulação 6-DOF de aeronave eVTOL tipo tailsitter biplace de 1500 kg em Python.**

Framework de simulação virtual orientado a objetos para veículos aéreos não-convencionais, com foco em **tailsitter** (4 rotores fixos com cant configurável) cobrindo todo o envelope: hover, transição, cruzeiro, retorno. Implementação direta da versão MATLAB original com revisão de engenharia.

![Status](https://img.shields.io/badge/status-functional-success)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## 📋 Descrição do Problema

### Aeronave alvo

Tailsitter biplace eVTOL com configuração não-convencional:
- **Massa**: 1500 kg (escala Joby S4 / Beta Alia)
- **Asa**: 11 m envergadura, 15 m² área, perfil NACA-0015 (simétrico bidirecional)
- **Propulsão**: 4 rotores fixos (raio 0.9 m, 4500 N cada) com **ângulos de cant configuráveis** em pitch e yaw
- **Empenagens**: superfícies de controle convencionais (elevator, aileron, rudder)
- **T/W**: 1.22 (margem para manobra)
- **V_LDmax**: 48.2 m/s a α=6.9°, L/D=15.6

### Missão típica

```
Altitude (m)
    1000  ═══════════════════
              |    \
              |     \  cruise @ V_LDmax
              |
       20  ═══
              |  hover 20s
        0  ──┘
              0          ~3 km North
```

### Desafios de engenharia

1. **Transição hover→cruzeiro**: rotação de 90° em pitch com controle de empuxo simultâneo
2. **Singularidade Euler em hover**: resolvido via quatérnios Hamilton
3. **Aerodinâmica 360°**: extrapolação Viterna-Corrigan para alto AoA
4. **Interação prop-asa**: slipstream sobre superfícies portantes
5. **Autoridade de atitude limitada**: ~8s mínimos para pitch 90°
6. **Cascade de controle robusto**: outer PD+FF + inner SO(3) geométrico

---

## 🏗 Arquitetura do Framework (5 Etapas)

| Etapa | Conteúdo |
|---|---|
| **1. Cinemática & Dinâmica** | Newton-Euler 6-DOF com quatérnios Hamilton + tensor inércia com $J_{xz}$ |
| **2. Aero & Propulsão** | Viterna-Corrigan 360° + BEMT + Glauert + slipstream + ESC 1ª ordem |
| **3. Controle** | Cascade: PD+FF outer (1 rad/s) + SO(3) Lee-Leok-McClamroch inner (8 rad/s) + alocador pseudoinversa |
| **4. Trajetória** | Planicidade diferencial + polinômios min-snap 8ª ordem + GA refinamento |
| **5. Ambiente** | RK4 500 Hz + Dryden MIL-F-8785C + ISA standard atmosphere |

---

## 🚀 Instalação

### Requisitos

- Python ≥ 3.10
- NumPy ≥ 1.24
- SciPy ≥ 1.10
- Matplotlib ≥ 3.7

### Instalação

```bash
git clone https://github.com/TobiasMorais/EVTOL-Tailsitter.git
cd EVTOL-Tailsitter
pip install -e .
```

Ou para desenvolvimento:

```bash
pip install -r requirements.txt
```

---

## 💻 Uso Rápido

### 1. Sanity Check Físico

Valida 23 indicadores físicos do modelo (T/W, L/D, V_stall, etc.):

```bash
python scripts/sanity_check.py
```

### 2. Suíte de Testes

15 suítes, 100+ testes unitários e de integração:

```bash
pytest tests/
```

### 3. Missão Completa

Roda missão realista (climb 20m → hover 20s → climb 1km → cruise) com plots:

```bash
python scripts/main.py
```

### 4. Uso em script Python

```python
from evtol.config import aircraft_config, controller_config, simulation_config
from evtol.simulation import run_simulation
from evtol.simulation.missions import build_realistic_mission
from evtol.visualization import plot_telemetry

# Setup
ac = aircraft_config()
ctrl = controller_config()
sim = simulation_config()

# Build mission
traj, info = build_realistic_mission(ac)
sim['t_final'] = traj.total_time() + 10

# Run
log = run_simulation(ac, ctrl, sim, traj)

# Visualize
plot_telemetry(log, ac)
```

---

## 📊 Resultados Principais

### Validação física (sanity_check)

| Métrica | Valor |
|---|---|
| T/W ratio | 1.22 ✓ |
| Hover throttle | 81.7% ✓ |
| Wing loading | 100 kg/m² ✓ |
| L/D max | 15.6 ✓ |
| V_LDmax | 48.2 m/s ✓ |
| Bandwidth cascade (outer/inner) | 8× ✓ |

### Testes (suíte completa)

```
[STAGE 1] Cinemática/Dinâmica       27/27 PASS
[STAGE 2] Aero-Propulsão            17/17 PASS
[STAGE 3] Controle                  33/33 PASS
[STAGE 4] Trajetória                15/15 PASS
[STAGE 5] Ambiente                  11/11 PASS
[STAGE 6] Integração (hover)         6/6 PASS  (drift 0.18 m em 10 s)
TOTAL: 109/109 PASS
```

---

## 📚 Documentação

| Documento | Conteúdo |
|---|---|
| [docs/TECHNICAL_REPORT.md](docs/TECHNICAL_REPORT.md) | Relatório técnico completo (13 seções, 21 refs) |
| [docs/EQUATIONS.md](docs/EQUATIONS.md) | Equações detalhadas das 5 etapas |
| [docs/MISSION_DESIGN.md](docs/MISSION_DESIGN.md) | Engenharia da missão (cálculos $V_{LDmax}$) |
| [docs/SANITY_CHECKS.md](docs/SANITY_CHECKS.md) | Validação física de parâmetros |

---

## 🧩 Estrutura do Repositório

```
EVTOL-Tailsitter/
├── README.md
├── pyproject.toml
├── requirements.txt
├── LICENSE
├── evtol/                          # pacote principal
│   ├── config/                     # parâmetros (aeronave, controle, simulação)
│   ├── core/                       # Aircraft, quaternion, SO(3)
│   ├── aerodynamics/               # Viterna-Corrigan + BEMT
│   ├── control/                    # SO(3) inner, PD+FF outer, allocator
│   ├── trajectory/                 # flatness, min-snap, GA
│   ├── environment/                # ISA, Dryden
│   ├── simulation/                 # RK4, simulator, missions
│   └── visualization/              # plots
├── scripts/
│   ├── main.py                     # missão completa
│   ├── sanity_check.py             # validação física
│   └── run_all_tests.py
├── tests/                          # ~15 suítes pytest
└── docs/                           # relatórios técnicos
```

---

## 🔬 Referências Acadêmicas

1. Lee, Leok, McClamroch (2010). *Geometric tracking control of a quadrotor UAV on SE(3)*. IEEE CDC.
2. Mellinger, Kumar (2011). *Minimum snap trajectory generation and control for quadrotors*. IEEE ICRA.
3. Viterna, Corrigan (1981). *Fixed Pitch Rotor Performance of Large Horizontal Axis Wind Turbines*. NASA CP-2230.
4. Leishman (2006). *Principles of Helicopter Aerodynamics*. Cambridge.
5. Smeur, Chu, de Croon (2016). *Adaptive Incremental Nonlinear Dynamic Inversion for Attitude Control of MAVs*. JGCD 39(3).
6. Stevens, Lewis, Johnson (2015). *Aircraft Control and Simulation*. Wiley.
7. MIL-F-8785C (1980). *Flying Qualities of Piloted Airplanes*.
8. Verling et al. (2016). *Full attitude control of a VTOL tailsitter UAV*. IEEE ICRA.

Lista completa em [docs/TECHNICAL_REPORT.md](docs/TECHNICAL_REPORT.md).

---

## 📄 Licença

MIT License — ver [LICENSE](LICENSE).

---

## 👤 Autor

**Tobias Morais** — Master's Research, eVTOL Control Systems

Originalmente desenvolvido em MATLAB, transcrito para Python com revisão de engenharia.

---

## 🤝 Contribuições

PRs são bem-vindos. Áreas de interesse:
- L1 adaptive augmentation para INDI
- VLM aerodynamics (substituir strip theory)
- Hardware-in-the-Loop (HIL) interface
- Modos de falha de rotor
- Sensores realistas + EKF

Issues: [github.com/TobiasMorais/EVTOL-Tailsitter/issues](https://github.com/TobiasMorais/EVTOL-Tailsitter/issues)
