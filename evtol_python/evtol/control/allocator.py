"""Control allocator: virtual moments -> actuator commands.

Maps:
    nu = [F_x_B; F_y_B; F_z_B; M_x; M_y; M_z]  -> u = [T_1, T_2, T_3, T_4, d_e, d_a, d_r]

via weighted pseudoinverse:
    u = W^-1 B^T (B W^-1 B^T)^-1 nu
"""
from __future__ import annotations

import numpy as np


class ControlAllocator:
    """Pseudoinverse allocator with saturation."""

    def __init__(self, ctrl_cfg: dict, ac_cfg: dict, rotor_axes: np.ndarray):
        self.cfg = ctrl_cfg
        self.ac_cfg = ac_cfg
        self.rotor_axes = rotor_axes
        self.B = self._build_effectiveness()

    def _build_effectiveness(self) -> np.ndarray:
        """Build 6 x 7 effectiveness matrix."""
        B = np.zeros((6, 7))
        positions = self.ac_cfg["rotor_position"]

        for i in range(self.ac_cfg["n_rotors"]):
            e_i = self.rotor_axes[:, i]
            r_i = positions[i]
            B[:3, i] = e_i
            B[3:, i] = np.cross(r_i, e_i)

        # Surface effectiveness (will be scaled by q_bar at allocation time)
        S = self.ac_cfg["htail"]["area"]
        arm = self.ac_cfg["htail"]["arm"]
        # Elevator: pitch moment (col 5, index 4)
        B[4, 4] = -self.ac_cfg["htail"]["elev_eff"] * S * arm
        # Aileron: roll moment (col 6, index 5)
        B[3, 5] = self.ac_cfg["ail"]["eff"] * self.ac_cfg["wing"]["area"] * self.ac_cfg["wing"]["span"]
        # Rudder: yaw moment (col 7, index 6)
        B[5, 6] = -self.ac_cfg["vtail"]["rud_eff"] * self.ac_cfg["vtail"]["area"] * self.ac_cfg["vtail"]["arm"]

        return B

    def allocate(self, nu: np.ndarray, q_bar: float) -> tuple[np.ndarray, dict]:
        """Allocate virtual command nu to actuators."""
        B_eff = self.B.copy()
        B_eff[:, 4:] *= max(q_bar, 1.0)  # scale surfaces by dynamic pressure

        W = self.cfg["alloc"]["W"]
        Wi = np.linalg.inv(W)
        BWB = B_eff @ Wi @ B_eff.T + 1e-3 * np.eye(6)
        u_unc = Wi @ B_eff.T @ np.linalg.solve(BWB, nu)

        # Saturate
        u = u_unc.copy()
        T_min = self.cfg["alloc"]["T_min"]
        T_max = self.cfg["alloc"]["T_max"]
        u[:4] = np.clip(u[:4], T_min, T_max)
        u[4] = np.clip(u[4], -self.ac_cfg["surf"]["elev_max"], self.ac_cfg["surf"]["elev_max"])
        u[5] = np.clip(u[5], -self.ac_cfg["ail"]["delta_max"], self.ac_cfg["ail"]["delta_max"])
        u[6] = np.clip(u[6], -self.ac_cfg["surf"]["rud_max"], self.ac_cfg["surf"]["rud_max"])

        info = {"unconstrained": u_unc, "residual": nu - B_eff @ u, "B": B_eff}
        return u, info

    def attitude_effectiveness(self) -> np.ndarray:
        """3 x N: d(omega_dot)/dT_i for rotors (used by INDI)."""
        n = self.ac_cfg["n_rotors"]
        G = np.zeros((3, n))
        positions = self.ac_cfg["rotor_position"]
        J_inv = self.ac_cfg["J_inv"]
        for i in range(n):
            M_per_T = np.cross(positions[i], self.rotor_axes[:, i])
            G[:, i] = J_inv @ M_per_T
        return G
