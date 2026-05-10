# Sanity Checks Físicos do Modelo eVTOL

Bateria de **18 verificações de plausibilidade física** que validam se o modelo se comporta consistentemente com primeiros princípios da aerodinâmica e mecânica de voo. Diferente dos testes unitários (que validam código), estes verificam se o **modelo físico** está coerente.

## Como rodar

```matlab
cd evtol_simulation
addpath(genpath(pwd))
results = sanity_check_model();
```

Saída: tabela com `[OK]`, `[WARN]`, `[FAIL]` e `[INFO]` para cada parâmetro derivado, com nota explicando a faixa esperada para a classe de aeronave.

## Lista das verificações

| # | Verificação | Fórmula | Faixa esperada | Razão física |
|---|---|---|---|---|
| 1 | Massa & peso | $m, W=mg$ | 1500 kg / 14715 N | Sanidade básica |
| 2 | T/W ratio máximo | $\sum T_{max} / W$ | 1.20 – 1.60 | Tailsitter precisa $>$ 1.20 (margem manobra) |
| 3 | Throttle de hover | $T_{hover}/T_{max}$ | 60% – 85% | Saudável: margem para controle |
| 4 | Wing loading | $W/S$ | 80–180 kg/m² | Faixa eVTOL leve |
| 5 | Disk loading | $W/(N\cdot A_{disk})$ | 100–500 kg/m² | Open rotor típico |
| 6 | Velocidade induzida hover | $v_i = \sqrt{T/(2\rho A)}$ | 10–30 m/s | Glauert |
| 7 | Potência de hover | $P = T\cdot v_i / FM$ | 100–500 kW | $FM \approx 0.7$ (Figure of Merit) |
| 8 | Power loading (kg/kW) | $m/P$ | 2–5 kg/kW | eVTOL saudável |
| 9 | Velocidade de estol | $V_s = \sqrt{2W/(\rho S C_{L,max})}$ | 25–40 m/s | Critério de design |
| 10 | Velocidade de cruzeiro | $1.3 V_s$ | 30–50 m/s | Margem de segurança |
| 11 | $C_L$ de cruzeiro | $2W/(\rho S V^2)$ | 0.3–0.8 | Abaixo de $C_{L,max}$ |
| 12 | $\alpha$ de cruzeiro | $C_L / C_{L_\alpha} + \alpha_0$ | 0–10° | Abaixo do estol (15°) |
| 13 | L/D de cruzeiro | $C_L/C_D$ | 8–15 | Aeronave fixa típica |
| 14 | Potência de cruzeiro | $D\cdot V$ | 30–150 kW | Razão hover/cruise 2-5× |
| 15 | $T = k_T \Omega^2$ | Estático | $\le 5\%$ erro | Calibração rotor |
| 16 | Separação de bandwidth | $\tau_{outer}/\tau_{inner}/\tau_{motor}$ | $\ge 1.5\times$ | Estabilidade cascata |
| 17 | $J_{yy}/J_{xx}$ | Inércia | 1.0–4.0 | Tailsitter com $J$ longitudinal maior |
| 18 | $J_{xz}/J_{yy}$ | Acoplamento | $< 0.20$ | Pequeno, não-dominante |
| 19 | Autoridade $M_y$ máxima | $\sum r_y \times T_{max}$ | $> J_{yy}\pi^2/2$ | Capaz de pitch 90°/1s |
| 20 | Slipstream far-wake | $V_{slip}=2 v_i$ | $\sim 2v_i$ | Teorema de Froude |
| 21 | Crossover aero/rotor | $V$ tal que $q\cdot S\cdot$arm$\cdot\delta = M_{rotor}$ | 10–30 m/s | Em cruzeiro, superfícies dominam |
| 22 | Queda livre $a = g$ | Sem F, M | $|a-g| < 10^{-9}$ | Conservação básica |
| 23 | Trim hover $T = mg$ | F_prop = mg | $|a| < 10^{-9}$ | Balanço estático |

## Interpretação dos resultados

### `[OK]`
Valor dentro da faixa típica esperada. Modelo coerente neste aspecto.

### `[WARN]`
Valor fora do esperado mas não absurdo. **Investigar** se intencional ou se algum parâmetro precisa de ajuste:
- Pode indicar configuração incomum (ex: aeronave overpowered ou underpowered)
- Pode indicar erro de unidade ou de fator (ex: $C_{L,max}$ subestimado)

### `[FAIL]`
Valor contradiz a física básica. **Bug no modelo**:
- T/W < 1: aeronave não consegue pairar
- $V_{stall} > V_{cruise}$: aeronave nunca atinge cruzeiro
- $L/D < 1$: drag maior que lift (impossível em cruzeiro estável)
- Queda livre com $|a-g| > 10^{-6}$: bug na dinâmica

### `[INFO]`
Valores diagnósticos sem critério pass/fail. Para inspeção visual.

## Como usar para validação de design

1. **Após cada alteração em `aircraft_config.m`**, rode `sanity_check_model()` para confirmar que os parâmetros derivados continuam coerentes.

2. **Comparar com aeronaves de referência**: Joby S4 (2200 kg, 5 rotores, ~600 kW hover), Beta Alia (2700 kg, 4 rotores), Lilium (3175 kg). Os valores derivados aqui devem estar na mesma ordem.

3. **Detecção de erro de unidade**: se algum `[WARN]` aparecer em vários parâmetros relacionados, provável erro de unidade (ex: chord em pés vs metros).

4. **Trade-off de design**: ajustar `aircraft_config.m` (ex: aumentar wing.area ou wing.CL_max) e ver impacto em $V_{stall}$, $L/D$, etc.
