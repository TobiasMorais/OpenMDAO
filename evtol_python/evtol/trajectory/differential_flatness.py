"""Differential flatness trajectory in (x, y, z, psi) with min-snap polynomials.

8th-order rest-to-rest polynomials per segment, C^3 continuous.
Solved in NORMALIZED time u = tau/T_seg to avoid numerical ill-conditioning
when T is large (T^7 = 10^14 for T=100s would give cond > 10^16).
"""
from __future__ import annotations

import numpy as np


class DifferentialFlatness:
    """Multi-segment min-snap polynomial trajectory."""

    def __init__(self, waypoints: np.ndarray, time_alloc: np.ndarray):
        """
        Parameters
        ----------
        waypoints : (N+1, 4) array of [x, y, z, psi] waypoints
        time_alloc : (N,) array of segment durations
        """
        self.waypoints = np.asarray(waypoints, dtype=float)
        self.time_alloc = np.asarray(time_alloc, dtype=float).flatten()
        self.coeffs = self._solve_minsnap()
        self.T_cumul = np.concatenate([[0.0], np.cumsum(self.time_alloc)])

    def _solve_minsnap(self) -> np.ndarray:
        """Solve min-snap polynomial system in NORMALIZED time u = tau/T.

        Returns coeffs of shape (4, 8, N_seg).
        Polynomial in u: p(u) = sum_{j=0..7} c_j * u^j
        Derivatives in physical time: dp/dtau = (dp/du)/T, etc.
        """
        N_seg = self.waypoints.shape[0] - 1
        C = np.zeros((4, 8, N_seg))

        # Normalized A matrix (entries in {0, 1, 2, 6, 12, 20, 30, 42, 60, 120, 210})
        A = np.array([
            [1, 0, 0, 0,  0,  0,   0,   0],  # p(0)
            [0, 1, 0, 0,  0,  0,   0,   0],  # p'(0)
            [0, 0, 2, 0,  0,  0,   0,   0],  # p''(0)
            [0, 0, 0, 6,  0,  0,   0,   0],  # p'''(0)
            [1, 1, 1, 1,  1,  1,   1,   1],  # p(1)
            [0, 1, 2, 3,  4,  5,   6,   7],  # p'(1)
            [0, 0, 2, 6, 12, 20,  30,  42],  # p''(1)
            [0, 0, 0, 6, 24, 60, 120, 210],  # p'''(1)
        ], dtype=float)

        for k in range(N_seg):
            p0 = self.waypoints[k]
            p1 = self.waypoints[k + 1]
            for d in range(4):
                b = np.array([p0[d], 0, 0, 0, p1[d], 0, 0, 0])
                C[d, :, k] = np.linalg.solve(A, b)
        return C

    def total_time(self) -> float:
        return float(self.T_cumul[-1])

    def eval(self, t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
        """Evaluate trajectory at time t. Returns (pos, vel, acc, jerk, snap, psi, psi_dot)."""
        t = float(np.clip(t, 0.0, self.T_cumul[-1]))

        # Find segment
        k = max(0, int(np.searchsorted(self.T_cumul, t, side="right")) - 1)
        k = min(k, self.coeffs.shape[2] - 1)
        T_seg = self.time_alloc[k]
        tau = t - self.T_cumul[k]
        u = tau / T_seg

        # Normalized basis vectors
        U = np.array([u ** j for j in range(8)])
        Ud = np.concatenate([[0.0], [j * u ** (j - 1) for j in range(1, 8)]])
        Udd = np.concatenate([[0.0, 0.0],
                              [j * (j - 1) * u ** (j - 2) for j in range(2, 8)]])
        Uddd = np.concatenate([[0.0, 0.0, 0.0],
                               [j * (j - 1) * (j - 2) * u ** (j - 3) for j in range(3, 8)]])
        Udddd = np.concatenate([[0.0, 0.0, 0.0, 0.0],
                                [j * (j - 1) * (j - 2) * (j - 3) * u ** (j - 4) for j in range(4, 8)]])

        sigma = np.zeros(4)
        sigmad = np.zeros(4)
        sigmadd = np.zeros(4)
        sigmaddd = np.zeros(4)
        sigmadddd = np.zeros(4)

        for d in range(4):
            c = self.coeffs[d, :, k]
            sigma[d] = c @ U
            sigmad[d] = (c @ Ud) / T_seg
            sigmadd[d] = (c @ Udd) / T_seg ** 2
            sigmaddd[d] = (c @ Uddd) / T_seg ** 3
            sigmadddd[d] = (c @ Udddd) / T_seg ** 4

        return sigma[:3], sigmad[:3], sigmadd[:3], sigmaddd[:3], sigmadddd[:3], float(sigma[3]), float(sigmad[3])
