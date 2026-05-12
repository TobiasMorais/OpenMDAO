"""Simulation configuration."""
from __future__ import annotations

import numpy as np


def simulation_config() -> dict:
    """Returns simulation parameters (solver, time, environment)."""
    sim: dict = {}

    sim["solver"]   = "rk4"
    sim["dt"]       = 0.002
    sim["t_final"]  = 600.0
    sim["scenario"] = "realistic"

    sim["env"] = {
        "gravity": 9.80665,
        "rho0":    1.225,
        "T0_K":    288.15,
        "h0":      0.0,
    }

    sim["wind"] = {
        "constant_NED": np.array([3.0, -1.0, 0.0]),
        "dryden": {
            "enable": True,
            "W20":    7.7,
            "h_ft":   100.0,
            "seed":   42,
        },
    }

    sim["sensors"] = {
        "imu": {
            "bias_gyro":    np.deg2rad([0.05, -0.03, 0.04]),
            "noise_gyro":   np.deg2rad(0.01),
            "bias_acc":     np.array([0.02, -0.01, 0.03]),
            "noise_acc":    0.05,
            "dt":           0.002,
            "enable_noise": True,
        },
    }

    return sim
