"""Stage 3: Control modules (SO3, PDFF, allocator, force_to_attitude)."""
import numpy as np

from evtol.config import aircraft_config, controller_config
from evtol.core import Aircraft, Quaternion
from evtol.core.so3 import SO3
from evtol.control.so3_controller import AttitudeControllerSO3
from evtol.control.pdff import PositionControllerPDFF
from evtol.control.allocator import ControlAllocator
from evtol.control.force_to_attitude import force_to_attitude


class TestSO3Controller:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        self.ctrl_cfg = controller_config()
        self.ctrl = AttitudeControllerSO3(self.ctrl_cfg)

    def test_equilibrium(self):
        M = self.ctrl.compute(np.eye(3), np.zeros(3), np.eye(3),
                              np.zeros(3), np.zeros(3), self.ac_cfg["J"], 0.0)
        assert np.linalg.norm(M) < 1e-9

    def test_bandwidth(self):
        for axis in range(3):
            omega_n = np.sqrt(self.ctrl_cfg["so3"]["kR"][axis, axis] / self.ac_cfg["J"][axis, axis])
            assert omega_n >= 5.0, f"axis {axis}: omega_n = {omega_n:.2f}"


class TestPDFF:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        self.ctrl_cfg = controller_config()
        self.pdff = PositionControllerPDFF(self.ctrl_cfg, self.ac_cfg["mass"])

    def test_hover_at_reference(self):
        p = np.array([0.0, 0.0, -10.0])
        v = np.zeros(3)
        p_ref = np.tile(p, (21, 1)).T
        v_ref = np.zeros((3, 21))
        f = self.pdff.compute(p, v, p_ref, v_ref)
        # At reference: f_cmd = -g_NED = [0, 0, -9.81]
        np.testing.assert_allclose(f, [0, 0, -9.80665], atol=1e-9)


class TestAllocator:
    def setup_method(self):
        self.ac_cfg = aircraft_config()
        self.ctrl_cfg = controller_config()
        ac = Aircraft(self.ac_cfg)
        self.alloc = ControlAllocator(self.ctrl_cfg, self.ac_cfg, ac.rotor_axes)

    def test_hover_share(self):
        W = self.ac_cfg["mass"] * 9.80665
        nu = np.array([W, 0, 0, 0, 0, 0])
        u, info = self.alloc.allocate(nu, 0.0)
        T_per_rotor = u[:4]
        assert np.all(T_per_rotor > 0)
        T_along_x = np.sum(T_per_rotor)
        assert abs(T_along_x - W) / W < 0.05


class TestForceToAttitude:
    def test_hover_x_B_up(self):
        F_des = np.array([0, 0, -9.80665])
        qd, T_total, Wd, Wd_dot = force_to_attitude(F_des, 0.0, 0.0, None)
        R = Quaternion.to_rotation_matrix(qd)
        x_B = R @ np.array([1, 0, 0])
        np.testing.assert_allclose(x_B, [0, 0, -1], atol=1e-6)

    def test_det_positive(self):
        qd, _, _, _ = force_to_attitude(np.array([0, 0.5, -9.81]), 0.0, 0.0, None)
        R = Quaternion.to_rotation_matrix(qd)
        assert abs(np.linalg.det(R) - 1) < 1e-6

    def test_heading_rotation(self):
        F = np.array([0, 0, -9.81])
        qd0, _, _, _ = force_to_attitude(F, 0.0, 0.0, None)
        qd1, _, _, _ = force_to_attitude(F, np.pi / 4, 0.0, None)
        y_B0 = Quaternion.to_rotation_matrix(qd0) @ np.array([0, 1, 0])
        y_B1 = Quaternion.to_rotation_matrix(qd1) @ np.array([0, 1, 0])
        assert np.linalg.norm(y_B0 - y_B1) > 0.1
