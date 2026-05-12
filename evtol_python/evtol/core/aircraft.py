"""Aircraft: 6-DOF rigid body dynamics with quaternion attitude.

State vector x (13 components):
    x[0:3]   = p_NED       position in inertial NED frame [m]
    x[3:6]   = v_NED       velocity in NED [m/s]
    x[6:10]  = q           quaternion body->world (Hamilton, scalar first)
    x[10:13] = omega_body  body angular rates [p; q; r] [rad/s]
"""
from __future__ import annotations

import numpy as np

from evtol.core.quaternion import Quaternion


class Aircraft:
    """Rigid-body 6-DOF tailsitter dynamics."""

    GRAVITY = 9.80665

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.rotor_axes = self.compute_rotor_axes()  # 3 x N body-frame thrust axes

    def compute_rotor_axes(self) -> np.ndarray:
        """Build unit thrust direction for each rotor (apply cant)."""
        n = self.cfg["n_rotors"]
        axes = np.zeros((3, n))
        for i in range(n):
            ay = self.cfg["rotor_yaw_cant"][i]
            ap = self.cfg["rotor_pitch_cant"][i]
            Rz = np.array([
                [np.cos(ay), -np.sin(ay), 0],
                [np.sin(ay),  np.cos(ay), 0],
                [0,           0,          1],
            ])
            Ry = np.array([
                [ np.cos(ap), 0, np.sin(ap)],
                [          0, 1,          0],
                [-np.sin(ap), 0, np.cos(ap)],
            ])
            axes[:, i] = Rz @ Ry @ np.array([1.0, 0.0, 0.0])
        return axes

    def dynamics(
        self,
        x: np.ndarray,
        F_aero_B: np.ndarray,
        M_aero_B: np.ndarray,
        F_prop_B: np.ndarray,
        M_prop_B: np.ndarray,
        M_gyro_B: np.ndarray,
    ) -> np.ndarray:
        """Newton-Euler 6-DOF dynamics. Returns x_dot (13,)."""
        v_NED = x[3:6]
        q = Quaternion.normalize(x[6:10])
        w_B = x[10:13]

        R_BW = Quaternion.to_rotation_matrix(q)  # body -> world

        # Translation (NED frame)
        F_total_B = F_aero_B + F_prop_B
        a_NED = (R_BW @ F_total_B) / self.cfg["mass"] + np.array([0, 0, self.GRAVITY])

        # Rotation (body frame), with J_xz coupling
        J = self.cfg["J"]
        M_total_B = M_aero_B + M_prop_B - M_gyro_B - np.cross(w_B, J @ w_B)
        wdot_B = self.cfg["J_inv"] @ M_total_B

        # Quaternion kinematics
        qdot = Quaternion.kinematic(q, w_B)

        return np.concatenate([v_NED, a_NED, qdot, wdot_B])

    def compute_rotor_gyro(self, omega_rotors_signed: np.ndarray, w_B: np.ndarray) -> np.ndarray:
        """M_gyro = omega_B x H_rotor."""
        J_r = self.cfg["rotor"]["J_rotor"]
        H = np.zeros(3)
        for i in range(self.cfg["n_rotors"]):
            H += J_r * omega_rotors_signed[i] * self.rotor_axes[:, i]
        return np.cross(w_B, H)

    def initial_state(self, mode: str = "hover") -> np.ndarray:
        """Build initial state for canonical scenarios."""
        if mode == "hover":
            p0 = np.array([0.0, 0.0, 0.0])
            q0 = Quaternion.from_euler(0, np.deg2rad(90), 0)
        elif mode == "hover_at_altitude":
            p0 = np.array([0.0, 0.0, -10.0])
            q0 = Quaternion.from_euler(0, np.deg2rad(90), 0)
        elif mode == "hover_at_altitude_init":
            p0 = np.array([0.0, 0.0, -1.0])
            q0 = Quaternion.from_euler(0, np.deg2rad(90), 0)
        elif mode == "cruise":
            p0 = np.array([0.0, 0.0, -50.0])
            q0 = Quaternion.from_euler(0, 0, 0)
            return np.concatenate([p0, np.array([50, 0, 0]), q0, np.zeros(3)])
        else:
            raise ValueError(f"Unknown init mode: {mode}")

        v0 = np.zeros(3)
        w0 = np.zeros(3)
        return np.concatenate([p0, v0, q0, w0])
