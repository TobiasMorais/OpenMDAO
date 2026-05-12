"""Geometric SO(3) attitude controller (Lee, Leok, McClamroch 2010)."""
from __future__ import annotations

import numpy as np

from evtol.core.so3 import SO3


class AttitudeControllerSO3:
    """Geometric attitude tracking on SO(3) with integrator and feed-forward.

    Control law (body frame):
        M = -k_R e_R - k_Omega e_Omega - k_I e_I
            + omega x J omega
            - J (hat(omega) R^T Rd Wd - R^T Rd Wd_dot)
    """

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.ei = np.zeros(3)
        self.prev_t: float | None = None

    def reset(self) -> None:
        self.ei[:] = 0
        self.prev_t = None

    def compute(
        self,
        R: np.ndarray,
        omega: np.ndarray,
        Rd: np.ndarray,
        Wd: np.ndarray,
        Wd_dot: np.ndarray,
        J: np.ndarray,
        t: float,
    ) -> np.ndarray:
        """Compute commanded body-frame moment."""
        e_R = SO3.err_mat(R, Rd)

        Wd_b = R.T @ Rd @ Wd
        e_W = omega - Wd_b

        # Integral with anti-windup
        if self.prev_t is None:
            dt = 0.0
        else:
            dt = max(0.0, t - self.prev_t)
        self.prev_t = t
        self.ei = self.ei + dt * e_R
        I_max = self.cfg["so3"]["I_max"]
        self.ei = np.clip(self.ei, -I_max, I_max)

        # Feed-forward (saturated to prevent FD-noise amplification)
        Wd_safe = np.clip(Wd, -3.0, 3.0)
        Wd_dot_safe = np.clip(Wd_dot, -2.0, 2.0)
        term_ff = -J @ (SO3.hat(omega) @ R.T @ Rd @ Wd_safe - R.T @ Rd @ Wd_dot_safe)

        M = (
            -self.cfg["so3"]["kR"] @ e_R
            - self.cfg["so3"]["kOm"] @ e_W
            - self.cfg["so3"]["kI"] @ self.ei
            + np.cross(omega, J @ omega)
            + term_ff
        )
        return M
