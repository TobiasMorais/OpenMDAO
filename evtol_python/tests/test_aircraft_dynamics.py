"""Stage 1c: 6-DOF rigid-body dynamics."""
import numpy as np

from evtol.config import aircraft_config
from evtol.core import Aircraft
from evtol.core.quaternion import Quaternion
from evtol.simulation.rk4 import rk4_step


class TestAircraftDynamics:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        self.ac = Aircraft(self.ac_cfg)

    def test_free_fall(self):
        x0 = np.concatenate([np.zeros(6), [1, 0, 0, 0], np.zeros(3)])
        xdot = self.ac.dynamics(x0, np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3))
        np.testing.assert_allclose(xdot[3:6], [0, 0, self.ac.GRAVITY], atol=1e-12)

    def test_quaternion_norm_preserved(self):
        x = self.ac.initial_state("hover")
        F_prop_B = np.array([self.ac_cfg["mass"] * self.ac.GRAVITY, 0, 0])
        odefun = lambda t, xx: self.ac.dynamics(xx, np.zeros(3), np.zeros(3), F_prop_B,
                                                  np.zeros(3), np.zeros(3))
        dt = 0.002
        for k in range(1000):
            x = rk4_step(odefun, k * dt, x, dt)
        assert abs(np.linalg.norm(x[6:10]) - 1) < 1e-6

    def test_hover_trim(self):
        x_h = self.ac.initial_state("hover")
        F_prop_B = np.array([self.ac_cfg["mass"] * self.ac.GRAVITY, 0, 0])
        xdot = self.ac.dynamics(x_h, np.zeros(3), np.zeros(3), F_prop_B, np.zeros(3), np.zeros(3))
        assert np.linalg.norm(xdot[3:6]) < 1e-9

    def test_pure_moment(self):
        x_h = self.ac.initial_state("hover")
        M_test = np.array([1.0, 0.0, 0.0])
        xdot = self.ac.dynamics(x_h, np.zeros(3), M_test, np.zeros(3), np.zeros(3), np.zeros(3))
        expected = self.ac_cfg["J_inv"] @ M_test
        np.testing.assert_allclose(xdot[10:13], expected, atol=1e-9)

    def test_rotor_gyro_paired_cancel(self):
        Omega_signed = np.array([200.0, -200.0, -200.0, 200.0])
        w_B = np.array([0.0, 1.0, 0.0])
        M = self.ac.compute_rotor_gyro(Omega_signed, w_B)
        assert np.linalg.norm(M) < 1e-3
