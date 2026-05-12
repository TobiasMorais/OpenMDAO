"""SO(3) Lie-algebra helpers."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


class SO3:
    """Static Lie-algebra operations on SO(3)."""

    @staticmethod
    def hat(v: ArrayLike) -> np.ndarray:
        """Skew-symmetric matrix [v]_x."""
        v = np.asarray(v, dtype=float).flatten()
        return np.array([
            [0,     -v[2],  v[1]],
            [v[2],   0,    -v[0]],
            [-v[1],  v[0],  0],
        ])

    @staticmethod
    def vee(X: ArrayLike) -> np.ndarray:
        """Inverse hat."""
        X = np.asarray(X, dtype=float)
        return np.array([X[2, 1], X[0, 2], X[1, 0]])

    @staticmethod
    def err_mat(R: ArrayLike, Rd: ArrayLike) -> np.ndarray:
        """Geometric attitude error: e_R = 0.5 vee(Rd^T R - R^T Rd).

        Lee, Leok, McClamroch (2010), Eq. (10).
        """
        R = np.asarray(R, dtype=float)
        Rd = np.asarray(Rd, dtype=float)
        S = 0.5 * (Rd.T @ R - R.T @ Rd)
        return SO3.vee(S)

    @staticmethod
    def exp(phi: ArrayLike) -> np.ndarray:
        """so(3) -> SO(3) (Rodrigues formula)."""
        phi = np.asarray(phi, dtype=float).flatten()
        a = np.linalg.norm(phi)
        if a < 1e-9:
            return np.eye(3) + SO3.hat(phi)
        ax = phi / a
        K = SO3.hat(ax)
        return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * (K @ K)

    @staticmethod
    def log(R: ArrayLike) -> np.ndarray:
        """SO(3) -> so(3)."""
        R = np.asarray(R, dtype=float)
        c = np.clip((np.trace(R) - 1) / 2, -1.0, 1.0)
        a = np.arccos(c)
        if a < 1e-9:
            return 0.5 * SO3.vee(R - R.T)
        return (a / (2 * np.sin(a))) * SO3.vee(R - R.T)
