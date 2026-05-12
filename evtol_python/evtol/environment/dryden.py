"""MIL-F-8785C Dryden turbulence model (low-altitude).

Linear filter formulation: drives unit-variance white noise through three
shaping filters (one per wind axis). Length scales and intensities follow
the standard low-altitude specification.
"""
from __future__ import annotations

import numpy as np


class DrydenWind:
    """Dryden turbulence + steady mean wind."""

    def __init__(self, cfg: dict):
        self.cfg = cfg
        seed = cfg["dryden"]["seed"]
        self.rng = np.random.default_rng(seed)
        self.x_u = 0.0
        self.x_v = np.zeros(2)
        self.x_w = np.zeros(2)
        self.constant_NED = np.asarray(cfg["constant_NED"], dtype=float)

    def step(self, V_airspeed: float, h_m: float, dt: float) -> np.ndarray:
        """Returns wind velocity in NED [m/s]."""
        if not self.cfg["dryden"]["enable"]:
            return self.constant_NED.copy()

        h_ft = max(50.0, h_m * 3.281)
        V = max(2.0, V_airspeed)
        W20 = self.cfg["dryden"]["W20"]

        denom = (0.177 + 0.000823 * h_ft) ** 1.2
        L_w = h_ft
        L_u = h_ft / denom
        L_v = L_u

        sigma_w = 0.1 * W20
        sigma_u = sigma_w / (0.177 + 0.000823 * h_ft) ** 0.4
        sigma_v = sigma_u

        # u channel (1st-order)
        tau_u = L_u / V
        a_u = np.exp(-dt / tau_u)
        n1 = self.rng.standard_normal()
        self.x_u = a_u * self.x_u + sigma_u * np.sqrt(1 - a_u ** 2) * n1
        u_g = self.x_u

        # v channel (2nd-order, simplified Tustin)
        tau_v = L_v / V
        a_v = np.exp(-dt / tau_v)
        n2 = self.rng.standard_normal()
        new_v1 = a_v * self.x_v[0] + sigma_v * np.sqrt(1 - a_v ** 2) * n2
        new_v2 = a_v * self.x_v[1] + (1 - a_v) * (
            new_v1 + np.sqrt(3) * tau_v * (new_v1 - self.x_v[0]) / dt
        )
        self.x_v = np.array([new_v1, new_v2])
        v_g = self.x_v[1]

        # w channel
        tau_w = L_w / V
        a_w = np.exp(-dt / tau_w)
        n3 = self.rng.standard_normal()
        new_w1 = a_w * self.x_w[0] + sigma_w * np.sqrt(1 - a_w ** 2) * n3
        new_w2 = a_w * self.x_w[1] + (1 - a_w) * (
            new_w1 + np.sqrt(3) * tau_w * (new_w1 - self.x_w[0]) / dt
        )
        self.x_w = np.array([new_w1, new_w2])
        w_g = self.x_w[1]

        return self.constant_NED + np.array([u_g, v_g, w_g])

    def reset(self) -> None:
        self.x_u = 0.0
        self.x_v[:] = 0
        self.x_w[:] = 0
