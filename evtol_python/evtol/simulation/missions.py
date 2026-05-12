"""Mission builders: realistic and hover scenarios."""
from __future__ import annotations

import numpy as np

from evtol.trajectory.differential_flatness import DifferentialFlatness


def build_hover_only(alt_m: float = 10.0, duration_s: float = 3.0) -> tuple[DifferentialFlatness, dict]:
    """Stationary hover trajectory (single segment, essentially constant)."""
    waypoints = np.array([
        [0.0, 0.0, -alt_m,        0.0],
        [0.0, 0.0, -alt_m - 0.001, 0.0],
    ])
    time_alloc = np.array([duration_s])
    traj = DifferentialFlatness(waypoints, time_alloc)
    info = {"init_mode": "hover_at_altitude", "alt_init_m": alt_m, "T_total": traj.total_time()}
    return traj, info


def build_realistic_mission(
    ac_cfg: dict,
    alt_init_m: float = 1.0,
    alt_low_m: float = 20.0,
    alt_cruise_m: float = 1000.0,
    cruise_distance_m: float = 1000.0,
    hover_duration_s: float = 20.0,
) -> tuple[DifferentialFlatness, dict]:
    """Build conservative 4-segment mission via DifferentialFlatness.

    Phases:
      1. Climb -1m -> -20m       (rest-to-rest)
      2. Hover at -20m for 20s   (zero-distance waypoint repeat)
      3. Climb -20m -> -1000m    (rest-to-rest, slow)
      4. Slow horizontal cruise  (rest-to-rest at 1000m)

    All phases at 30% of thrust envelope with 2x safety factor.
    The aircraft remains in tailsitter orientation throughout climbs.
    """
    g = 9.80665
    W = ac_cfg["mass"] * g
    rho = 1.225

    # Max L/D reference (informational)
    AR = ac_cfg["wing"]["AR"]
    e = ac_cfg["wing"]["e"]
    CD0 = ac_cfg["wing"]["CD0"]
    K = 1 / (np.pi * AR * e)
    CL_LDmax = np.sqrt(CD0 / K)
    V_LDmax = np.sqrt(2 * W / (rho * ac_cfg["wing"]["area"] * CL_LDmax))
    alpha_LDmax = CL_LDmax / ac_cfg["wing"]["CL_alpha"] + ac_cfg["wing"]["alpha0"]
    LD_max = 0.5 / np.sqrt(CD0 * K)

    # Conservative envelope (30%, 2x safety)
    a_envelope = (ac_cfg["n_rotors"] * ac_cfg["rotor"]["thrust_max"]) / ac_cfg["mass"]
    a_vert_max = max(0.3 * (a_envelope - g), 0.2)
    a_horiz_max = max(0.3 * np.sqrt(max(a_envelope ** 2 - g ** 2, 1.0)), 0.5)
    coef = 15

    T_climb1 = max(8.0, 2 * np.sqrt(coef * (alt_low_m - alt_init_m) / a_vert_max))
    T_hover = hover_duration_s
    d_alt = alt_cruise_m - alt_low_m
    T_climb2 = max(45.0, 2 * np.sqrt(coef * d_alt / a_vert_max))
    T_cruise = max(30.0, 2 * np.sqrt(coef * cruise_distance_m / a_horiz_max))

    waypoints = np.array([
        [0.0,               0.0, -alt_init_m,   0.0],   # start (above ground)
        [0.0,               0.0, -alt_low_m,    0.0],   # after climb to 20m
        [0.0,               0.0, -alt_low_m,    0.0],   # after hover
        [0.0,               0.0, -alt_cruise_m, 0.0],   # after climb to 1000m
        [cruise_distance_m, 0.0, -alt_cruise_m, 0.0],   # after cruise
    ])
    time_alloc = np.array([T_climb1, T_hover, T_climb2, T_cruise])
    traj = DifferentialFlatness(waypoints, time_alloc)

    info = {
        "V_LDmax": V_LDmax,
        "alpha_LDmax_deg": np.rad2deg(alpha_LDmax),
        "LD_max": LD_max,
        "CL_LDmax": CL_LDmax,
        "T_climb1": T_climb1,
        "T_hover": T_hover,
        "T_climb2": T_climb2,
        "T_cruise": T_cruise,
        "T_total": traj.total_time(),
        "distance_north_total": cruise_distance_m,
        "init_mode": "hover_at_altitude_init",
        "alt_init_m": alt_init_m,
    }

    print(f"""
========================================================
   Realistic Mission (DifferentialFlatness engine)
========================================================
   Aerodynamic design point (informational):
     V_LDmax = {V_LDmax:.1f} m/s ({V_LDmax*3.6:.0f} km/h)
     alpha_LDmax = {np.rad2deg(alpha_LDmax):.1f} deg
     L/D_max = {LD_max:.1f}

   Conservative trajectory (30% envelope, 2x safety):
     a_vert_max = {a_vert_max:.2f} m/s^2
     a_horiz_max = {a_horiz_max:.2f} m/s^2

   Initial state: aircraft at altitude {alt_init_m:.0f} m
   Phase durations:
     1. Climb {alt_init_m:.0f}->{alt_low_m:.0f}m    : {T_climb1:.1f} s
     2. Hover at {alt_low_m:.0f} m       : {T_hover:.1f} s
     3. Climb {alt_low_m:.0f}->{alt_cruise_m:.0f}m  : {T_climb2:.1f} s
     4. Cruise {cruise_distance_m:.0f}m N     : {T_cruise:.1f} s
     TOTAL                  : {info['T_total']:.1f} s ({info['T_total']/60:.2f} min)
========================================================
""")
    return traj, info
