"""Stage 2a: Viterna-Corrigan aerodynamics."""
import numpy as np

from evtol.config import aircraft_config
from evtol.aerodynamics import Aerodynamics


class TestAerodynamics:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        self.aero = Aerodynamics(self.ac_cfg)
        self.defs = self.aero.build_surface_strips()

    def test_zero_airspeed_zero_force(self):
        N = len(self.defs)
        F, M = self.aero.compute_forces_moments(
            np.zeros(3), np.zeros(3), self.defs, np.zeros((3, N)), 1.225
        )
        assert np.linalg.norm(F) + np.linalg.norm(M) < 1e-6

    def test_symmetric_airfoil_zero_alpha(self):
        CL, CD = self.aero.airfoil_360(0)
        assert abs(CL) < 0.05
        assert abs(CD - self.ac_cfg["wing"]["CD0"]) < 0.005

    def test_linear_regime(self):
        for alpha_deg in [-5, -2, 0, 2, 5, 8]:
            alpha = np.deg2rad(alpha_deg)
            CL, _ = self.aero.airfoil_360(alpha)
            expected = self.ac_cfg["wing"]["CL_alpha"] * alpha
            assert abs(CL - expected) < 0.05

    def test_post_stall_drops(self):
        CL15, _ = self.aero.airfoil_360(np.deg2rad(15))
        CL20, _ = self.aero.airfoil_360(np.deg2rad(20))
        assert abs(CL20) < abs(CL15)

    def test_flat_plate_90deg(self):
        _, CD = self.aero.airfoil_360(np.deg2rad(90))
        assert 0.8 < abs(CD) < 1.5

    def test_symmetric_lift(self):
        CLp, _ = self.aero.airfoil_360(np.deg2rad(8))
        CLn, _ = self.aero.airfoil_360(np.deg2rad(-8))
        assert abs(CLp + CLn) < 0.05

    def test_dynamic_pressure_scaling(self):
        N = len(self.defs)
        F1, _ = self.aero.compute_forces_moments(
            np.array([20.0, 0, 0]), np.zeros(3), self.defs, np.zeros((3, N)), 1.225
        )
        F2, _ = self.aero.compute_forces_moments(
            np.array([40.0, 0, 0]), np.zeros(3), self.defs, np.zeros((3, N)), 1.225
        )
        ratio = np.linalg.norm(F2) / max(np.linalg.norm(F1), 1e-9)
        assert abs(ratio - 4) < 0.5

    def test_slipstream_increases_force(self):
        N = len(self.defs)
        V_slip = np.zeros((3, N))
        V_slip[:, 0] = [25, 0, 0]
        F_no, _ = self.aero.compute_forces_moments(
            np.zeros(3), np.zeros(3), self.defs, np.zeros((3, N)), 1.225
        )
        F_yes, _ = self.aero.compute_forces_moments(
            np.zeros(3), np.zeros(3), self.defs, V_slip, 1.225
        )
        assert np.linalg.norm(F_yes) > np.linalg.norm(F_no)
