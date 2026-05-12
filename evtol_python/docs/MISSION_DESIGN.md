# Realistic Mission Design — Engineering Rationale

## Mission Profile

```
   Altitude (m)
        |
   1000 |═════════════════════════════════════════════
        |       /
        |      /  (Phase 3: climb-and-transition)
        |     /
        |    /
        |   /
        |  /
        | /
     20 |═══════ (Phase 2: hover 20s)
        |/
        | (Phase 1: vertical climb)
      0 +----------------------------> North (m)
        0                            ~3 km
```

## Phase Specifications

### Phase 1 — Vertical Climb (0 → 20 m)
- **Type**: rest-to-rest 8th-order min-snap polynomial
- **Initial state**: at origin, hovering, $\theta = +90°$ (tailsitter upright)
- **Final state**: at $-20$ m altitude, zero velocity
- **Duration**: ~14 s (limited by 50% of vertical thrust envelope, $a_{max,vert} \approx 1.1$ m/s²)
- **Pitch profile**: maintains $90°$ throughout (pure vertical thrust)
- **Tracking expectation**: error $< 5$ m (mostly initial transient from motor spinup)

### Phase 2 — Hover Hold (20 m, 20 s)
- **Type**: stationary
- **Position**: $[0, 0, -20]$ m, fixed
- **Duration**: 20 s (per user spec)
- **Aircraft state**: $\theta = +90°$, all velocities zero
- **Tracking expectation**: error $< 1$ m (cascade has 5 s settling time)
- **Purpose**: validate steady-state stability

### Phase 3 — Climb-and-Transition (20 → 1000 m)
- **Type**: rest-to-cruise 8th-order polynomial (asymmetric BCs)
- **Initial state**: at $-20$ m, hovering
- **Final state**: at $-1000$ m altitude AND moving at $V_{LDmax}$ northward
- **Duration**: ~120 s (limited by climb rate)
- **Pitch profile**: smoothly transitions from $90°$ to $\alpha_{LDmax} \approx 7°$
- **Distance traveled**: ~2.9 km north (during the climb-acceleration)
- **Engineering challenge**: most demanding phase — body must rotate $83°$ while accelerating to $V_{LDmax}$ and gaining $980$ m of altitude
- **Tracking expectation**: error grows during transition, then settles in cruise

### Phase 4 — Steady Cruise (1000 m, ~3 km)
- **Type**: constant velocity
- **Velocity**: $[V_{LDmax}, 0, 0]$ = $[48.3, 0, 0]$ m/s
- **Duration**: $3000 / 48.3 \approx 62$ s
- **Pitch**: $\alpha_{LDmax} \approx 7°$ (max L/D angle of attack)
- **Tracking expectation**: error $< 30$ m steady state

## Mathematical Design — Max L/D Cruise Point

For a parabolic drag polar:
$$C_D = C_{D_0} + K C_L^2, \quad K = \frac{1}{\pi AR \cdot e}$$

Maximum $L/D$ occurs when:
$$\frac{d(L/D)}{d\alpha} = 0 \;\Rightarrow\; C_{D_0} = K C_L^2$$
$$\boxed{C_{L,LDmax} = \sqrt{C_{D_0}/K}}$$

Required cruise velocity from lift balance:
$$L = W \;\Rightarrow\; \tfrac{1}{2}\rho V^2 S \cdot C_{L,LDmax} = W$$
$$\boxed{V_{LDmax} = \sqrt{\frac{2W}{\rho S \sqrt{C_{D_0}/K}}}}$$

Maximum efficiency:
$$\boxed{\left(\frac{L}{D}\right)_{max} = \frac{1}{2\sqrt{C_{D_0} K}}}$$

For our 1500 kg tailsitter:
- $C_{D_0} = 0.022, AR = 8.07, e = 0.85$
- $K = 0.0467$
- $C_{L,LDmax} = 0.687$
- $V_{LDmax} = 48.3$ m/s ($174$ km/h)
- $\alpha_{LDmax} = 6.9°$
- $L/D_{max} = 15.6$

## MissionTrajectory Class Architecture

```matlab
traj = MissionTrajectory({...
    MissionTrajectory.make_rest_to_rest(p_origin, p_alt_low, T1),
    MissionTrajectory.make_hover(p_alt_low, T_hover),
    MissionTrajectory.make_rest_to_cruise(p_alt_low, p_cruise_start, v_cruise, T_transition),
    MissionTrajectory.make_cruise(p_cruise_start, v_cruise, T_cruise)
});
```

Each phase has:
- `type`: one of `'hover'`, `'cruise'`, `'rest_to_rest'`, `'rest_to_cruise'`, `'cruise_to_rest'`
- `duration`: time span [s]
- Position/velocity boundary conditions
- For polynomial phases: 8 coefficients per axis (8th order satisfies 4 BCs at each end)

The `eval(t)` method dispatches to the correct phase, returning $(p, v, a, j, s, \psi, \dot\psi)$ continuous up to the third derivative.

## Validation Test (Stage 7)

`test_integration_realistic.m` runs the full mission and validates per-phase metrics:

| Check | Threshold | Justification |
|---|---|---|
| Phase 1 (climb) max error | < 20 m | Transient from initial state |
| Phase 2 (hover) max error | < 5 m | Steady-state hover |
| Phase 3 (transition) max error | < 200 m | Most demanding; allows transient |
| Phase 4 (cruise settled) max error | < 100 m | After 30 s settling |
| Cruise pitch | within 30° of $\alpha_{LDmax}$ | Aerodynamic trim |
| Quaternion norm | $< 10^{-3}$ error | Numerical integrity |
| Final altitude | within 200 m of 1000 m | No catastrophic divergence |

## Engineering Decisions

1. **Why 50% thrust envelope?** Reserves 50% control margin so NMPC and SO(3) have authority to compensate disturbances and unmodeled effects.

2. **Why $V_{LDmax}$ as cruise speed?** Maximum aerodynamic efficiency = minimum power = maximum range/endurance per unit energy.

3. **Why min-snap 8th-order polynomial?** $C^3$ continuous (snap = 4th derivative) trajectories produce smooth thrust commands without high-frequency content that would saturate actuators.

4. **Why two separate stages 1 and 3?** The hover-to-cruise transition is operationally and physically distinct from pure vertical climb. Splitting allows independent tuning of climb rate and transition aggression.

5. **Why $20$ s hover hold?** Per user requirement. Also serves as cascade convergence test before complex maneuvers.

## Future Extensions

- Phase 5: cruise-to-rest deceleration to landing-prep hover
- Phase 6: vertical descent and touchdown
- Wind disturbance during cruise (Dryden enabled)
- Energy-optimal phase 3 (currently max-margin; GA could optimize)
- Multi-waypoint cruise (zigzag, climb-and-glide)
