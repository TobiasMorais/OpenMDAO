"""Aircraft configuration: 1500 kg biplace tailsitter eVTOL."""
from __future__ import annotations

import numpy as np

from typing import Any


def aircraft_config() -> dict[str, Any]:
    """Returns aircraft parameter dict (Joby S4/Beta Alia-scaled tailsitter)."""
    ac: dict[str, Any] = {}

    # --- Mass & inertia ---
    ac["mass"] = 1500.0
    ac["cg_body"] = np.array([0.05, 0.0, 0.0])

    # Inertia tensor (kg*m^2), plane of symmetry X-Z
    ac["J"] = np.array([
        [2200.0,    0.0,  180.0],
        [   0.0, 3800.0,    0.0],
        [ 180.0,    0.0, 5500.0],
    ])
    ac["J_inv"] = np.linalg.inv(ac["J"])

    # --- Wing ---
    ac["wing"] = {
        "span":         11.0,
        "chord":         1.4,
        "area":         15.0,
        "AR":           11.0 ** 2 / 15.0,    # 8.07
        "e":             0.85,
        "airfoil":      "NACA-0015",
        "alpha0":        0.0,
        "CL_alpha":      5.7,
        "CL_max":        1.35,
        "alpha_stall":   np.deg2rad(15.0),
        "CD0":           0.022,
        "Cm0":           0.0,
        "Cm_alpha":     -0.45,
    }

    # --- Horizontal tail ---
    ac["htail"] = {
        "area":      2.5,
        "arm":       4.5,
        "CL_alpha":  4.5,
        "eta":       0.9,
        "elev_eff":  0.55,
    }

    # --- Vertical tail ---
    ac["vtail"] = {
        "area":     1.6,
        "arm":      4.5,
        "CY_beta": -0.55,
        "rud_eff":  0.45,
    }

    # --- Control surfaces ---
    ac["ail"] = {"eff": 0.18, "delta_max": np.deg2rad(25.0)}
    ac["surf"] = {
        "elev_max": np.deg2rad(25.0),
        "rud_max":  np.deg2rad(25.0),
        "tau":      0.04,
    }

    # --- Rotors (4 fixed, 2 per wingtip) ---
    ac["n_rotors"] = 4

    # Position vectors (body frame: X fwd, Y right, Z down)
    ac["rotor_position"] = np.array([
        [0.30, +5.0, -0.15],   # rotor 1: right wingtip, upper
        [0.30, +5.0, +0.15],   # rotor 2: right wingtip, lower
        [0.30, -5.0, -0.15],   # rotor 3: left  wingtip, upper
        [0.30, -5.0, +0.15],   # rotor 4: left  wingtip, lower
    ])

    # Rotation direction (+1 CCW viewed from front, -1 CW)
    ac["rotor_spin_dir"] = np.array([+1, -1, -1, +1])

    # Cant angles (parametric sweep)
    ac["rotor_pitch_cant"] = np.deg2rad([+2.0, -2.0, +2.0, -2.0])
    ac["rotor_yaw_cant"]   = np.deg2rad([+1.0, +1.0, -1.0, -1.0])

    # Rotor performance
    rotor = {
        "radius":          0.90,
        "n_blades":        5,
        "blade_chord":     0.10,
        "blade_pitch":     np.deg2rad(12.0),
        "cl_alpha_blade":  5.7,
        "cd0_blade":       0.012,
        "J_rotor":         0.32,
        "omega_max":     230.0,
        "thrust_max":   4500.0,
        "tau_motor":       0.06,
    }
    rotor["disk_area"] = np.pi * rotor["radius"] ** 2
    rotor["solidity"] = rotor["n_blades"] * rotor["blade_chord"] / (np.pi * rotor["radius"])
    rotor["kT"] = rotor["thrust_max"] / rotor["omega_max"] ** 2
    rotor["kQ"] = 0.025 * rotor["kT"] * rotor["radius"]
    ac["rotor"] = rotor

    # Slipstream geometry
    ac["slipstream"] = {
        "bath_matrix": np.array([
            [1, 0, 1, 0],   # rotor 1 -> [wing_R, wing_L, htail, vtail]
            [1, 0, 1, 1],
            [0, 1, 1, 0],
            [0, 1, 1, 1],
        ]),
        "k_s":   1.0,
        "x_off": 0.8,
    }

    # Battery
    ac["battery"] = {"capacity_kWh": 280.0, "mass_kg": 420.0}

    return ac
