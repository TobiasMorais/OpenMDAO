"""Stage 2b: Propulsion (BEMT + ESC + slipstream)."""
import numpy as np

from evtol.config import aircraft_config
from evtol.core import Aircraft
from evtol.aerodynamics import Propulsion


class TestPropulsion:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        ac = Aircraft(self.ac_cfg)
        self.prop = Propulsion(self.ac_cfg, ac.rotor_axes)

    def test_zero_speed_zero_force(self):
        F, M, _, _ = self.prop.compute(np.zeros(3), 1.225)
        assert np.linalg.norm(F) + np.linalg.norm(M) < 1e-9

    def test_T_kT_omega2(self):
        self.prop.Omega_actual = 0.5 * self.ac_cfg["rotor"]["omega_max"] * np.ones(4)
        F, _, _, _ = self.prop.compute(np.zeros(3), 1.225)
        T_expected = 4 * self.ac_cfg["rotor"]["kT"] * (0.5 * self.ac_cfg["rotor"]["omega_max"]) ** 2
        T_along_x = F[0]
        assert abs(T_along_x - T_expected) / T_expected < 0.02

    def test_esc_convergence(self):
        from evtol.core import Aircraft
        ac = Aircraft(self.ac_cfg)
        prop = Propulsion(self.ac_cfg, ac.rotor_axes)
        T_cmd = 2000.0 * np.ones(4)
        dt = 0.001
        N = int(5 * self.ac_cfg["rotor"]["tau_motor"] / dt)
        for _ in range(N):
            prop.step_actuator(T_cmd, dt)
        T_actual = self.ac_cfg["rotor"]["kT"] * prop.Omega_actual ** 2
        assert np.all(T_actual > 0.98 * T_cmd)

    def test_saturation(self):
        prop = Propulsion(self.ac_cfg, np.eye(3, 4))
        T_huge = 1e5 * np.ones(4)
        for _ in range(100):
            prop.step_actuator(T_huge, 0.001)
        T_sat = self.ac_cfg["rotor"]["kT"] * prop.Omega_actual ** 2
        assert np.all(T_sat <= self.ac_cfg["rotor"]["thrust_max"] + 1)

    def test_hover_induced_velocity(self):
        from evtol.core import Aircraft
        ac = Aircraft(self.ac_cfg)
        prop = Propulsion(self.ac_cfg, ac.rotor_axes)
        T_hover = self.ac_cfg["mass"] * 9.80665 / 4
        prop.Omega_actual = np.sqrt(T_hover / self.ac_cfg["rotor"]["kT"]) * np.ones(4)
        _, _, vi, _ = prop.compute(np.zeros(3), 1.225)
        v_expected = np.sqrt(T_hover / (2 * 1.225 * self.ac_cfg["rotor"]["disk_area"]))
        assert abs(vi[0] - v_expected) / v_expected < 0.02

    def test_forward_flight_decreases_vi(self):
        from evtol.core import Aircraft
        ac = Aircraft(self.ac_cfg)
        prop = Propulsion(self.ac_cfg, ac.rotor_axes)
        T_hover = self.ac_cfg["mass"] * 9.80665 / 4
        prop.Omega_actual = np.sqrt(T_hover / self.ac_cfg["rotor"]["kT"]) * np.ones(4)
        _, _, vi_hover, _ = prop.compute(np.zeros(3), 1.225)
        _, _, vi_fwd, _ = prop.compute(np.array([30, 0, 0]), 1.225)
        assert vi_fwd[0] < vi_hover[0]

    def test_hover_trim_equality(self):
        from evtol.core import Aircraft
        ac = Aircraft(self.ac_cfg)
        prop = Propulsion(self.ac_cfg, ac.rotor_axes)
        T_per = self.ac_cfg["mass"] * 9.80665 / 4
        prop.Omega_actual = np.sqrt(T_per / self.ac_cfg["rotor"]["kT"]) * np.ones(4)
        F, _, _, _ = prop.compute(np.zeros(3), 1.225)
        W = self.ac_cfg["mass"] * 9.80665
        assert abs(F[0] - W) / W < 0.05
