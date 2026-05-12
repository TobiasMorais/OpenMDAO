"""PD with acceleration feedforward (industry-standard outer loop).

Control law (NED frame):
    f_cmd = -K_p (p - p_ref) - K_v (v - v_ref) + a_ref - g_NED

Reference: Mellinger & Kumar (2011), Bouabdallah (2007), Lee-Leok-McClamroch (2010).
"""
from __future__ import annotations

import numpy as np


class PositionControllerPDFF:
    """Stateless PD + acceleration feedforward outer loop."""

    def __init__(self, cfg: dict, mass: float):
        self.cfg = cfg
        # Tuning: omega_n = 1 rad/s, zeta = 0.85
        omega_n = 1.0
        zeta = 0.85
        self.K_p = omega_n ** 2 * np.eye(3)        # 1.0 * I
        self.K_v = 2 * zeta * omega_n * np.eye(3)  # 1.7 * I
        self.a_max = cfg["nmpc"]["acc_max"]

    def reset(self) -> None:
        pass

    def compute(
        self,
        p: np.ndarray,
        v: np.ndarray,
        p_ref_traj: np.ndarray,  # 3 x (N+1)
        v_ref_traj: np.ndarray,
        a_ref: np.ndarray | None = None,
    ) -> np.ndarray:
        """Compute commanded specific force (NED, m/s^2)."""
        if a_ref is None:
            a_ref = np.zeros(3)

        p_ref = p_ref_traj[:, 0]
        v_ref = v_ref_traj[:, 0]

        err_p = p - p_ref
        err_v = v - v_ref
        g_NED = np.array([0.0, 0.0, 9.80665])

        f_cmd = -self.K_p @ err_p - self.K_v @ err_v + a_ref - g_NED

        # Saturate at a_max (direction preserved)
        mag = np.linalg.norm(f_cmd)
        if mag > self.a_max:
            f_cmd = f_cmd * (self.a_max / mag)
        return f_cmd
