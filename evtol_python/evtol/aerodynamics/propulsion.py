"""Propulsion: BEMT-derived induced velocity + 1st-order ESC dynamics.

For each rotor i:
    Actuator (ESC):  tau_m * dOmega_i/dt + Omega_i = Omega_cmd_i
    Static thrust:   T_i = k_T * Omega_i^2
    Static torque:   Q_i = k_Q * Omega_i^2
    Induced velocity (Glauert axial):
        v_i = -V_x/2 + sqrt((V_x/2)^2 + T/(2*rho*A))
"""
from __future__ import annotations

import numpy as np


class Propulsion:
    """4-rotor propulsion with ESC dynamics, BEMT, and slipstream."""

    def __init__(self, cfg: dict, rotor_axes: np.ndarray):
        self.cfg = cfg
        self.rotor_axes = rotor_axes  # 3 x N
        self.Omega_actual = np.zeros(cfg["n_rotors"])

    def step_actuator(self, T_cmd: np.ndarray, dt: float) -> None:
        """Advance ESC dynamics by dt."""
        T_cmd = np.clip(T_cmd, 0, self.cfg["rotor"]["thrust_max"])
        Omega_cmd = np.sqrt(T_cmd / self.cfg["rotor"]["kT"])
        tau = self.cfg["rotor"]["tau_motor"]
        self.Omega_actual = self.Omega_actual + (dt / tau) * (Omega_cmd - self.Omega_actual)
        self.Omega_actual = np.clip(self.Omega_actual, 0, self.cfg["rotor"]["omega_max"])

    def compute(self, V_B: np.ndarray, rho: float) -> tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray
    ]:
        """Returns (F_prop_B, M_prop_B, v_i_per_rotor, Omega_signed)."""
        n = self.cfg["n_rotors"]
        F_prop_B = np.zeros(3)
        M_prop_B = np.zeros(3)
        v_i = np.zeros(n)
        Omega_signed = np.zeros(n)

        kT = self.cfg["rotor"]["kT"]
        kQ = self.cfg["rotor"]["kQ"]
        A = self.cfg["rotor"]["disk_area"]
        spin_dir = self.cfg["rotor_spin_dir"]
        positions = self.cfg["rotor_position"]

        for i in range(n):
            Om = self.Omega_actual[i]
            T = kT * Om ** 2
            Q = kQ * Om ** 2

            axis_i = self.rotor_axes[:, i]
            r_i = positions[i]

            V_x = V_B @ axis_i

            if T > 1e-3:
                v_i[i] = -V_x / 2 + np.sqrt((V_x / 2) ** 2 + T / (2 * rho * A))
            Omega_signed[i] = spin_dir[i] * Om

            F_i = T * axis_i
            M_thrust = np.cross(r_i, F_i)
            M_reaction = -spin_dir[i] * Q * axis_i

            F_prop_B += F_i
            M_prop_B += M_thrust + M_reaction

        return F_prop_B, M_prop_B, v_i, Omega_signed

    def slipstream_velocities(self, v_i: np.ndarray, V_B: np.ndarray) -> np.ndarray:
        """Distributes induced velocity into 3 x N_surfaces array."""
        n_rotors = self.cfg["n_rotors"]
        bath = self.cfg["slipstream"]["bath_matrix"]
        n_surf = bath.shape[1]
        V_slip = np.zeros((3, n_surf))
        k_s = self.cfg["slipstream"]["k_s"]

        for s in range(n_surf):
            v_axial = 0.0
            contributors = 0
            for i in range(n_rotors):
                if bath[i, s] > 0:
                    v_axial += k_s * v_i[i]
                    contributors += 1
            if contributors > 0:
                V_slip[0, s] = v_axial / contributors
        return V_slip
