"""Stages 4 & 5: Trajectory, RK4, Dryden, atmosphere."""
import numpy as np

from evtol.trajectory import DifferentialFlatness
from evtol.simulation.rk4 import rk4_step
from evtol.environment import atmosphere_isa, DrydenWind


class TestDifferentialFlatness:
    def test_endpoints(self):
        wp = np.array([[0, 0, 0, 0], [10, 5, -20, 0], [30, 10, -50, 0]])
        times = np.array([3.0, 4.0])
        df = DifferentialFlatness(wp, times)
        p0, *_ = df.eval(0)
        np.testing.assert_allclose(p0, wp[0, :3], atol=1e-6)
        p_end, *_ = df.eval(sum(times))
        np.testing.assert_allclose(p_end, wp[2, :3], atol=1e-6)
        p_mid, *_ = df.eval(times[0])
        np.testing.assert_allclose(p_mid, wp[1, :3], atol=1e-6)

    def test_rest_to_rest_velocities(self):
        wp = np.array([[0, 0, 0, 0], [10, 5, -20, 0]])
        df = DifferentialFlatness(wp, np.array([3.0]))
        _, v0, *_ = df.eval(0)
        _, vT, *_ = df.eval(3.0)
        assert np.linalg.norm(v0) < 1e-6
        assert np.linalg.norm(vT) < 1e-6

    def test_zero_accel_at_endpoints(self):
        wp = np.array([[0, 0, 0, 0], [10, 5, -20, 0]])
        df = DifferentialFlatness(wp, np.array([3.0]))
        _, _, a0, *_ = df.eval(0)
        _, _, aT, *_ = df.eval(3.0)
        assert np.linalg.norm(a0) < 1e-6
        assert np.linalg.norm(aT) < 1e-6


class TestRK4:
    def test_exp_decay(self):
        """dx/dt = -x, x(0) = 1, RK4 should match exp(-t) to 4th order."""
        x = 1.0
        dt = 0.01
        N = 100
        for k in range(N):
            x = rk4_step(lambda t, xx: -xx, k * dt, x, dt)
        assert abs(x - np.exp(-1.0)) < 1e-8

    def test_harmonic_oscillator(self):
        """Energy conservation over one period."""
        state = np.array([1.0, 0.0])
        ode = lambda t, s: np.array([s[1], -s[0]])
        dt = 0.001
        N = int(2 * np.pi / dt)
        E0 = 0.5 * (state[0] ** 2 + state[1] ** 2)
        for k in range(N):
            state = rk4_step(ode, k * dt, state, dt)
        E1 = 0.5 * (state[0] ** 2 + state[1] ** 2)
        assert abs(E1 - E0) / E0 < 1e-3


class TestAtmosphere:
    def test_sea_level(self):
        assert abs(atmosphere_isa(0) - 1.225) < 1e-3

    def test_tropopause(self):
        rho_top = atmosphere_isa(11000)
        assert 0.3 < rho_top < 0.45

    def test_monotonic(self):
        heights = [0, 500, 1500, 3000, 6000, 10000]
        rhos = [atmosphere_isa(h) for h in heights]
        diffs = np.diff(rhos)
        assert np.all(diffs < 0)


class TestDryden:
    def test_disabled_returns_constant(self):
        cfg = {
            "constant_NED": np.array([3, -1, 0]),
            "dryden": {"enable": False, "W20": 7.7, "seed": 42},
        }
        wm = DrydenWind(cfg)
        v = wm.step(20.0, 100.0, 0.01)
        np.testing.assert_allclose(v, cfg["constant_NED"], atol=1e-15)

    def test_variance(self):
        cfg = {
            "constant_NED": np.array([0.0, 0, 0]),
            "dryden": {"enable": True, "W20": 7.7, "seed": 42},
        }
        wm = DrydenWind(cfg)
        N = 5000
        samples = np.array([wm.step(25.0, 50.0, 0.01) for _ in range(N)]).T
        sigma = np.std(samples, axis=1)
        assert np.all(sigma > 0.1)
