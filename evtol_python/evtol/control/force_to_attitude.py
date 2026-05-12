"""Force vector to desired attitude (SE(3) consistent for tailsitter).

Construction:
    x_B = F_des / |F_des|                 body x = thrust direction
    c_psi = [cos(psi); sin(psi); 0]       heading reference vector (NED ground plane)
    y_B = (c_psi x x_B) / |c_psi x x_B|   body y, perpendicular to x_B and c_psi
    z_B = x_B x y_B                       body z completes right-handed frame
"""
from __future__ import annotations

import numpy as np

from evtol.core.quaternion import Quaternion


def force_to_attitude(
    F_des_NED: np.ndarray,
    psi_d: float,
    dt: float,
    prev: dict | None,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """Map F_des to (qd, T_total, Wd, Wd_dot).

    Wd, Wd_dot estimated by finite differences between successive calls.
    """
    F_des_NED = np.asarray(F_des_NED).flatten()
    Fmag = np.linalg.norm(F_des_NED)

    if Fmag < 1e-3:
        qd = Quaternion.from_euler(0, np.deg2rad(90), psi_d)
        return qd, 0.0, np.zeros(3), np.zeros(3)

    x_B = F_des_NED / Fmag
    c_psi = np.array([np.cos(psi_d), np.sin(psi_d), 0.0])

    # y_B = (c_psi x x_B) / norm
    y_B_unnorm = np.cross(c_psi, x_B)
    if np.linalg.norm(y_B_unnorm) < 1e-3:
        # Singular: thrust parallel to heading vector. Fallback.
        if abs(x_B[2]) < 0.99:
            y_B_unnorm = np.cross(np.array([0.0, 0.0, 1.0]), x_B)
        else:
            y_B_unnorm = np.cross(np.array([0.0, 1.0, 0.0]), x_B)
    y_B = y_B_unnorm / np.linalg.norm(y_B_unnorm)
    z_B = np.cross(x_B, y_B)
    z_B = z_B / np.linalg.norm(z_B)

    R = np.column_stack([x_B, y_B, z_B])
    if np.linalg.det(R) < 0:
        z_B = -z_B
        R = np.column_stack([x_B, y_B, z_B])

    qd = Quaternion.from_rotation_matrix(R)
    T_total = Fmag

    # Wd, Wd_dot via finite difference
    if prev is None or prev.get("qd") is None or dt <= 0:
        return qd, T_total, np.zeros(3), np.zeros(3)

    qrel = Quaternion.mul(Quaternion.conj(prev["qd"]), qd)
    if qrel[0] < 0:
        qrel = -qrel
    phi = Quaternion.log_map(qrel)
    Wd = phi / dt

    if prev.get("Wd") is None:
        Wd_dot = np.zeros(3)
    else:
        Wd_dot = (Wd - prev["Wd"]) / dt
        Wd_dot = np.clip(Wd_dot, -50, 50)

    return qd, T_total, Wd, Wd_dot
