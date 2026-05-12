"""Stage 1a: Hamilton quaternion unit tests."""
import numpy as np
import pytest

from evtol.core.quaternion import Quaternion


class TestQuaternion:
    def test_identity_mul(self):
        q = Quaternion.normalize([0.5, 0.3, -0.7, 0.2])
        np.testing.assert_allclose(Quaternion.mul(q, Quaternion.identity()), q, atol=1e-12)

    def test_conjugate_inverse(self):
        q = Quaternion.normalize([0.5, 0.3, -0.7, 0.2])
        qe = Quaternion.mul(q, Quaternion.conj(q))
        np.testing.assert_allclose(qe, Quaternion.identity(), atol=1e-12)

    def test_R_identity(self):
        np.testing.assert_allclose(Quaternion.to_rotation_matrix(Quaternion.identity()),
                                   np.eye(3), atol=1e-12)

    def test_R_orthogonal(self):
        q = Quaternion.normalize([0.5, 0.3, -0.7, 0.2])
        R = Quaternion.to_rotation_matrix(q)
        np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-10)

    def test_euler_roundtrip(self):
        rng = np.random.default_rng(42)
        for _ in range(20):
            phi = -np.pi + 2 * np.pi * rng.random()
            th = np.deg2rad(-85 + 170 * rng.random())
            psi = -np.pi + 2 * np.pi * rng.random()
            q = Quaternion.from_euler(phi, th, psi)
            eu = Quaternion.to_euler(q)
            np.testing.assert_allclose(eu, [phi, th, psi], atol=1e-9)

    def test_kinematic_zero(self):
        q = Quaternion.normalize([0.5, 0.3, -0.7, 0.2])
        qd = Quaternion.kinematic(q, np.zeros(3))
        np.testing.assert_allclose(qd, np.zeros(4), atol=1e-15)

    def test_kinematic_z_axis(self):
        qd = Quaternion.kinematic(Quaternion.identity(), [0, 0, 1])
        np.testing.assert_allclose(qd, [0, 0, 0, 0.5], atol=1e-12)

    def test_exp_zero(self):
        np.testing.assert_allclose(Quaternion.exp_map(np.zeros(3)),
                                   Quaternion.identity(), atol=1e-15)

    def test_log_exp_inverse(self):
        phi = np.array([0.3, -0.5, 0.2])
        np.testing.assert_allclose(Quaternion.log_map(Quaternion.exp_map(phi)), phi, atol=1e-9)

    def test_rotation_z90(self):
        q = Quaternion.from_euler(0, 0, np.pi / 2)
        v = Quaternion.rotate(q, [1, 0, 0])
        np.testing.assert_allclose(v, [0, 1, 0], atol=1e-12)

    def test_err_mul_identity(self):
        q = Quaternion.normalize([0.5, 0.3, -0.7, 0.2])
        qe = Quaternion.err_mul(q, q)
        np.testing.assert_allclose(qe, Quaternion.identity(), atol=1e-12)

    def test_norm_preservation(self):
        rng = np.random.default_rng(7)
        q = Quaternion.identity()
        for _ in range(100):
            qr = Quaternion.normalize(rng.standard_normal(4))
            q = Quaternion.normalize(Quaternion.mul(q, qr))
        assert abs(np.linalg.norm(q) - 1) < 1e-9
