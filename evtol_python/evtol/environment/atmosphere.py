"""ISA standard atmosphere (troposphere)."""
from __future__ import annotations

import numpy as np


def atmosphere_isa(h_m: float, sim_env: dict | None = None) -> float:
    """Density via ISA standard atmosphere (troposphere only).

    T(h) = T0 - L * h    (L = 0.0065 K/m)
    p(h) = p0 * (T/T0)^(g/(R*L))
    rho  = p / (R * T)
    """
    if sim_env is None:
        T0 = 288.15
        rho0 = 1.225
    else:
        T0 = sim_env.get("T0_K", 288.15)
        rho0 = sim_env.get("rho0", 1.225)

    L = 0.0065
    g = 9.80665
    R = 287.058
    h = float(np.clip(h_m, 0.0, 11000.0))
    T = T0 - L * h
    return rho0 * (T / T0) ** (g / (R * L) - 1)
