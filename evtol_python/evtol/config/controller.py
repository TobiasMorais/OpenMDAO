"""Controller configuration."""
from __future__ import annotations

import numpy as np


def controller_config() -> dict:
    """Returns control gains and parameters."""
    ctrl: dict = {}

    # --- Inner attitude loop ---
    ctrl["inner"] = {"type": "SO3"}  # 'SO3' | 'INDI'

    # Geometric SO(3) gains scaled by inertia
    # omega_n = 8 rad/s, zeta = 0.7
    omega_n = 8.0
    zeta_att = 0.7
    J_diag = np.array([2200.0, 3800.0, 5500.0])
    ctrl["so3"] = {
        "kR":   np.diag(J_diag * omega_n ** 2),
        "kOm":  np.diag(2 * zeta_att * omega_n * J_diag),
        "kI":   np.diag(0.05 * J_diag),
        "I_max": np.deg2rad(15.0),
    }

    # INDI parameters (Smeur 2016)
    ctrl["indi"] = {
        "G":         None,                       # filled at init
        "lp_omega":  50.0,
        "lp_u":      50.0,
        "K_rate":    np.diag([20, 25, 12]),
        "K_att":     np.diag([10, 12,  6]),
    }

    # --- Outer position/velocity loop ---
    ctrl["outer"] = {"type": "PDFF"}  # 'PDFF' (default) | 'NMPC'

    # NMPC parameters (optional)
    ctrl["nmpc"] = {
        "N":         20,
        "dt":        0.05,
        "Q_pos":     np.diag([10, 10, 25]),
        "Q_vel":     np.diag([2, 2, 4]),
        "R_acc":     0.05 * np.eye(3),
        "Qf":        50 * np.eye(6),
        "acc_max":   2.5 * 9.81,
        "tilt_max":  np.deg2rad(85),
    }

    # --- Control allocator ---
    ctrl["alloc"] = {
        "method": "wpinv",
        "W":      np.diag([1, 1, 1, 1, 5, 5, 5]),
        "T_min":  0.0,
        "T_max":  4500.0,
        "dT_max": 50000.0,
    }

    # --- Loop frequencies ---
    ctrl["f_inner"] = 500.0
    ctrl["f_outer"] = 50.0

    return ctrl
