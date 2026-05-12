"""Hamilton quaternion utilities (scalar-first, body-to-world convention).

Singularity-free attitude representation. Operations follow the convention:
    q = [q0, q1, q2, q3]^T
    R(q) = (q0^2 - qv^T qv) I + 2 qv qv^T - 2 q0 [qv]_x  (body -> world)
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


class Quaternion:
    """Static Hamilton quaternion operations."""

    @staticmethod
    def identity() -> np.ndarray:
        return np.array([1.0, 0.0, 0.0, 0.0])

    @staticmethod
    def normalize(q: ArrayLike) -> np.ndarray:
        q = np.asarray(q, dtype=float).flatten()
        n = np.linalg.norm(q)
        if n < 1e-12:
            return Quaternion.identity()
        return q / n

    @staticmethod
    def conj(q: ArrayLike) -> np.ndarray:
        q = np.asarray(q, dtype=float).flatten()
        return np.array([q[0], -q[1], -q[2], -q[3]])

    @staticmethod
    def mul(a: ArrayLike, b: ArrayLike) -> np.ndarray:
        """Hamilton product a (x) b."""
        a = np.asarray(a, dtype=float).flatten()
        b = np.asarray(b, dtype=float).flatten()
        a0, av = a[0], a[1:]
        b0, bv = b[0], b[1:]
        return np.concatenate([
            [a0 * b0 - av @ bv],
            a0 * bv + b0 * av + np.cross(av, bv),
        ])

    @staticmethod
    def to_rotation_matrix(q: ArrayLike) -> np.ndarray:
        """Quaternion -> 3x3 rotation matrix (body -> world)."""
        q0, qx, qy, qz = np.asarray(q, dtype=float).flatten()
        return np.array([
            [1 - 2 * (qy**2 + qz**2),
             2 * (qx * qy - qz * q0),
             2 * (qx * qz + qy * q0)],
            [2 * (qx * qy + qz * q0),
             1 - 2 * (qx**2 + qz**2),
             2 * (qy * qz - qx * q0)],
            [2 * (qx * qz - qy * q0),
             2 * (qy * qz + qx * q0),
             1 - 2 * (qx**2 + qy**2)],
        ])

    @staticmethod
    def from_rotation_matrix(R: ArrayLike) -> np.ndarray:
        """3x3 rotation matrix -> quaternion."""
        R = np.asarray(R, dtype=float)
        tr = np.trace(R)
        if tr > 0:
            S = 2 * np.sqrt(tr + 1.0)
            q0 = 0.25 * S
            qx = (R[2, 1] - R[1, 2]) / S
            qy = (R[0, 2] - R[2, 0]) / S
            qz = (R[1, 0] - R[0, 1]) / S
        elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
            S = 2 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
            q0 = (R[2, 1] - R[1, 2]) / S
            qx = 0.25 * S
            qy = (R[0, 1] + R[1, 0]) / S
            qz = (R[0, 2] + R[2, 0]) / S
        elif R[1, 1] > R[2, 2]:
            S = 2 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
            q0 = (R[0, 2] - R[2, 0]) / S
            qx = (R[0, 1] + R[1, 0]) / S
            qy = 0.25 * S
            qz = (R[1, 2] + R[2, 1]) / S
        else:
            S = 2 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
            q0 = (R[1, 0] - R[0, 1]) / S
            qx = (R[0, 2] + R[2, 0]) / S
            qy = (R[1, 2] + R[2, 1]) / S
            qz = 0.25 * S
        return Quaternion.normalize([q0, qx, qy, qz])

    @staticmethod
    def from_euler(phi: float, theta: float, psi: float) -> np.ndarray:
        """Z-Y-X Tait-Bryan Euler -> quaternion."""
        cphi, sphi = np.cos(phi / 2), np.sin(phi / 2)
        cth, sth = np.cos(theta / 2), np.sin(theta / 2)
        cpsi, spsi = np.cos(psi / 2), np.sin(psi / 2)
        return np.array([
            cphi * cth * cpsi + sphi * sth * spsi,
            sphi * cth * cpsi - cphi * sth * spsi,
            cphi * sth * cpsi + sphi * cth * spsi,
            cphi * cth * spsi - sphi * sth * cpsi,
        ])

    @staticmethod
    def to_euler(q: ArrayLike) -> np.ndarray:
        """Quaternion -> Z-Y-X Tait-Bryan Euler angles [phi, theta, psi]."""
        q0, qx, qy, qz = np.asarray(q, dtype=float).flatten()
        phi = np.arctan2(2 * (q0 * qx + qy * qz), 1 - 2 * (qx**2 + qy**2))
        sinp = np.clip(2 * (q0 * qy - qz * qx), -1.0, 1.0)
        theta = np.arcsin(sinp)
        psi = np.arctan2(2 * (q0 * qz + qx * qy), 1 - 2 * (qy**2 + qz**2))
        return np.array([phi, theta, psi])

    @staticmethod
    def kinematic(q: ArrayLike, omega_body: ArrayLike) -> np.ndarray:
        """q_dot = 0.5 * Omega(omega_body) * q."""
        q = np.asarray(q, dtype=float).flatten()
        w = np.asarray(omega_body, dtype=float).flatten()
        Omega = np.array([
            [0,    -w[0], -w[1], -w[2]],
            [w[0],  0,     w[2], -w[1]],
            [w[1], -w[2],  0,     w[0]],
            [w[2],  w[1], -w[0],  0],
        ])
        return 0.5 * Omega @ q

    @staticmethod
    def err_mul(q: ArrayLike, qd: ArrayLike) -> np.ndarray:
        """Multiplicative error: qe = qd^-1 (x) q, with unwinding fix."""
        qe = Quaternion.mul(Quaternion.conj(qd), q)
        if qe[0] < 0:
            qe = -qe
        return qe

    @staticmethod
    def rotate(q: ArrayLike, v: ArrayLike) -> np.ndarray:
        """Rotate vector v from body to world by quaternion q."""
        return Quaternion.to_rotation_matrix(q) @ np.asarray(v, dtype=float).flatten()

    @staticmethod
    def exp_map(phi: ArrayLike) -> np.ndarray:
        """so(3) -> SU(2): exp(phi/2)."""
        phi = np.asarray(phi, dtype=float).flatten()
        a = np.linalg.norm(phi)
        if a < 1e-9:
            return np.concatenate([[1.0], 0.5 * phi])
        return np.concatenate([[np.cos(a / 2)], np.sin(a / 2) / a * phi])

    @staticmethod
    def log_map(q: ArrayLike) -> np.ndarray:
        """SU(2) -> so(3): 2*log(q). Returns rotation vector phi."""
        q = np.asarray(q, dtype=float).flatten()
        if q[0] < 0:
            q = -q
        v = q[1:]
        s = q[0]
        nv = np.linalg.norm(v)
        if nv < 1e-9:
            return 2 * v
        return 2 * np.arctan2(nv, s) * (v / nv)
