"""Stage 1b: SO(3) Lie-algebra unit tests."""
import numpy as np

from evtol.core.so3 import SO3


class TestSO3:
    def test_vee_hat_inverse(self):
        v = np.array([0.7, -0.2, 0.5])
        np.testing.assert_allclose(SO3.vee(SO3.hat(v)), v, atol=1e-15)

    def test_hat_antisymmetric(self):
        v = np.array([0.7, -0.2, 0.5])
        H = SO3.hat(v)
        np.testing.assert_allclose(H + H.T, np.zeros((3, 3)), atol=1e-15)

    def test_exp_zero(self):
        np.testing.assert_allclose(SO3.exp(np.zeros(3)), np.eye(3), atol=1e-15)

    def test_exp_log_inverse(self):
        rng = np.random.default_rng(11)
        for _ in range(5):
            phi = rng.standard_normal(3) * 0.7
            R = SO3.exp(phi)
            R_back = SO3.exp(SO3.log(R))
            np.testing.assert_allclose(R, R_back, atol=1e-9)

    def test_err_mat_same_R(self):
        R = SO3.exp(np.array([0.3, 0.1, -0.2]))
        e = SO3.err_mat(R, R)
        assert np.linalg.norm(e) < 1e-12

    def test_det_exp_is_one(self):
        phi = np.array([0.4, -0.3, 0.6])
        R = SO3.exp(phi)
        assert abs(np.linalg.det(R) - 1) < 1e-12

    def test_rodrigues_axis_z(self):
        ang = np.pi / 3
        R = SO3.exp(np.array([0, 0, 1]) * ang)
        Rexpected = np.array([
            [np.cos(ang), -np.sin(ang), 0],
            [np.sin(ang),  np.cos(ang), 0],
            [0, 0, 1],
        ])
        np.testing.assert_allclose(R, Rexpected, atol=1e-12)
