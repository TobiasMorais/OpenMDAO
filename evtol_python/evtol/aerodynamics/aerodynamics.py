"""Aerodynamics: strip-theory with Viterna-Corrigan 360-deg airfoil model.

Decomposes aircraft into 4 lifting surfaces (right wing, left wing, htail, vtail)
and applies post-stall extrapolation per Viterna-Corrigan (1981):
    C_L = A1 sin(2 alpha) + A2 cos^2(alpha) / sin(alpha)
    C_D = B1 sin^2(alpha) + B2 cos(alpha)
"""
from __future__ import annotations

import numpy as np


class Aerodynamics:
    """Strip-theory aerodynamics with 360-deg AoA support."""

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.viterna = self._precompute_viterna()

    def _precompute_viterna(self) -> dict:
        """Compute Viterna A1, A2, B1, B2 coefficients."""
        AR = self.cfg["wing"]["AR"]
        CD_max = 1.11 + 0.018 * AR  # Flat plate at 90 deg
        as_ = self.cfg["wing"]["alpha_stall"]
        CL_s = self.cfg["wing"]["CL_alpha"] * (as_ - self.cfg["wing"]["alpha0"])
        CD_s = self.cfg["wing"]["CD0"] + (CL_s ** 2) / (np.pi * AR * self.cfg["wing"]["e"])

        A1 = CD_max / 2
        A2 = (CL_s - CD_max * np.sin(as_) * np.cos(as_)) * np.sin(as_) / np.cos(as_) ** 2
        B1 = CD_max
        B2 = (CD_s - CD_max * np.sin(as_) ** 2) / np.cos(as_)

        return {"A1": A1, "A2": A2, "B1": B1, "B2": B2,
                "CD_max": CD_max, "alpha_stall": as_}

    def airfoil_360(self, alpha: float) -> tuple[float, float]:
        """Piecewise CL, CD with smooth transitions over alpha in [-pi, pi]."""
        alpha = np.arctan2(np.sin(alpha), np.cos(alpha))  # wrap
        as_ = self.viterna["alpha_stall"]

        if abs(alpha) <= as_:
            CL = self.cfg["wing"]["CL_alpha"] * (alpha - self.cfg["wing"]["alpha0"])
            CD = self.cfg["wing"]["CD0"] + (CL ** 2) / (
                np.pi * self.cfg["wing"]["AR"] * self.cfg["wing"]["e"]
            )
        else:
            a_abs = abs(alpha)
            if a_abs > np.pi / 2:
                a_abs = np.pi - a_abs  # mirror for back-side flow
            CL_v = self.viterna["A1"] * np.sin(2 * a_abs) + \
                   self.viterna["A2"] * np.cos(a_abs) ** 2 / max(np.sin(a_abs), 1e-3)
            CD_v = self.viterna["B1"] * np.sin(a_abs) ** 2 + self.viterna["B2"] * np.cos(a_abs)
            CL = np.sign(alpha) * CL_v
            CD = max(0.01, CD_v)

        # Hard saturation
        CL = float(np.clip(CL, -2.0, 2.0))
        CD = float(np.clip(CD, 0.0, 2.5))
        return CL, CD

    def build_surface_strips(self) -> list[dict]:
        """Build strip definitions for compute_forces_moments."""
        ac = self.cfg
        return [
            {"position": np.array([0.05, +ac["wing"]["span"] / 4, 0.0]),
             "area": ac["wing"]["area"] / 2, "role": "wing_R"},
            {"position": np.array([0.05, -ac["wing"]["span"] / 4, 0.0]),
             "area": ac["wing"]["area"] / 2, "role": "wing_L"},
            {"position": np.array([-ac["htail"]["arm"], 0.0, 0.0]),
             "area": ac["htail"]["area"], "role": "htail"},
            {"position": np.array([-ac["vtail"]["arm"], 0.0, -0.5]),
             "area": ac["vtail"]["area"], "role": "vtail"},
        ]

    def compute_forces_moments(
        self,
        V_B: np.ndarray,
        w_B: np.ndarray,
        surf_def: list[dict],
        V_slip_per_surf: np.ndarray,  # 3 x N_surfaces
        rho: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute aerodynamic force and moment in body frame."""
        F_B = np.zeros(3)
        M_B = np.zeros(3)

        for k, strip in enumerate(surf_def):
            r = strip["position"]
            S = strip["area"]
            role = strip["role"]

            V_local = V_B + np.cross(w_B, r) + V_slip_per_surf[:, k]
            V_mag = np.linalg.norm(V_local)
            if V_mag < 0.5:
                continue

            u, v_y, w = V_local
            alpha = np.arctan2(-w, u)
            beta = np.arcsin(np.clip(v_y / V_mag, -1.0, 1.0))

            CL, CD = self.airfoil_360(alpha)
            qbar = 0.5 * rho * V_mag ** 2
            L = qbar * S * CL
            D = qbar * S * CD

            if role in ("wing_R", "wing_L", "htail"):
                e_drag = -V_local / V_mag
                Vxz = np.array([V_local[0], 0.0, V_local[2]])
                nVxz = max(np.linalg.norm(Vxz), 1e-6)
                e_lift = np.array([-Vxz[2], 0.0, Vxz[0]]) / nVxz
                F_strip = L * e_lift + D * e_drag
            elif role == "vtail":
                Y = qbar * S * self.cfg["vtail"]["CY_beta"] * beta
                F_strip = np.array([-D, Y, 0.0])
            else:
                continue

            F_B += F_strip
            M_B += np.cross(r, F_strip)

        # Static pitching moment about CG
        V_inf_mag = np.linalg.norm(V_B)
        if V_inf_mag > 0.5:
            qbar = 0.5 * rho * V_inf_mag ** 2
            u, _, w = V_B
            alpha = np.arctan2(-w, u)
            Cm = self.cfg["wing"]["Cm0"] + self.cfg["wing"]["Cm_alpha"] * alpha
            M_B[1] += qbar * self.cfg["wing"]["area"] * self.cfg["wing"]["chord"] * Cm

        return F_B, M_B
