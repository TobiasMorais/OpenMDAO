"""Classical fixed-step Runge-Kutta 4 integrator with quaternion renormalization."""
from __future__ import annotations

import numpy as np
from typing import Callable


def rk4_step(
    odefun: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    x: np.ndarray,
    dt: float,
) -> np.ndarray:
    """RK4 step with quaternion renormalization (assumes q at x[6:10])."""
    k1 = odefun(t, x)
    k2 = odefun(t + dt / 2, x + dt / 2 * k1)
    k3 = odefun(t + dt / 2, x + dt / 2 * k2)
    k4 = odefun(t + dt, x + dt * k3)
    x_next = x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    if len(x_next) >= 10:
        q = x_next[6:10]
        nq = np.linalg.norm(q)
        if nq > 1e-9:
            x_next[6:10] = q / nq
        else:
            x_next[6:10] = np.array([1.0, 0.0, 0.0, 0.0])
    return x_next
