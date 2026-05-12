function tr = test_flatness()
% TEST_FLATNESS  Stage 4a — Differential flatness min-snap polynomials.
tr = test_helpers('init', 'test_flatness');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'trajectory'));

% Two-segment trajectory
wp = [ 0,  0,    0, 0;
      10,  5,  -20, 0;
      30, 10,  -50, 0];
times = [3.0; 4.0];
df = DifferentialFlatness(wp, times);

% T1: endpoints match
[p0, ~, ~, ~, ~, ~, ~] = df.eval(0.0);
[p_end, ~, ~, ~, ~, ~, ~] = df.eval(sum(times));
[p_mid, ~, ~, ~, ~, ~, ~] = df.eval(times(1));
tr = test_helpers('assert_near', tr, p0, wp(1,1:3)', 1e-6, 'T1a start point');
tr = test_helpers('assert_near', tr, p_end, wp(3,1:3)', 1e-6, 'T1b end point');
tr = test_helpers('assert_near', tr, p_mid, wp(2,1:3)', 1e-6, 'T1c mid point');

% T2: rest-to-rest velocities zero
[~, v0, ~, ~, ~, ~, ~] = df.eval(0.0);
[~, vT, ~, ~, ~, ~, ~] = df.eval(sum(times));
tr = test_helpers('assert_lt', tr, norm(v0), 1e-6, 'T2a v(0) = 0');
tr = test_helpers('assert_lt', tr, norm(vT), 1e-6, 'T2b v(T) = 0');

% T3: rest-to-rest accelerations zero
[~, ~, a0, ~, ~, ~, ~] = df.eval(0.0);
[~, ~, aT, ~, ~, ~, ~] = df.eval(sum(times));
tr = test_helpers('assert_lt', tr, norm(a0), 1e-6, 'T3a a(0) = 0');
tr = test_helpers('assert_lt', tr, norm(aT), 1e-6, 'T3b a(T) = 0');

% T4: rest-to-rest jerk zero
[~, ~, ~, j0, ~, ~, ~] = df.eval(0.0);
[~, ~, ~, jT, ~, ~, ~] = df.eval(sum(times));
tr = test_helpers('assert_lt', tr, norm(j0), 1e-6, 'T4a jerk(0) = 0');
tr = test_helpers('assert_lt', tr, norm(jT), 1e-6, 'T4b jerk(T) = 0');

% T5: continuity at segment boundary (position + velocity)
eps_t = 1e-6;
[p_minus, v_minus, ~, ~, ~, ~, ~] = df.eval(times(1) - eps_t);
[p_plus,  v_plus,  ~, ~, ~, ~, ~] = df.eval(times(1) + eps_t);
tr = test_helpers('assert_near', tr, p_minus, p_plus, 1e-3, 'T5a position continuous at junction');
tr = test_helpers('assert_near', tr, v_minus, v_plus, 1.0, 'T5b velocity continuous at junction (loose due to FD)');

% T6: thrust direction recovery — at mid-segment, x_B should align with m*(acc-g)
[~, ~, acc, ~, ~, ~, ~] = df.eval(0.5 * sum(times));
F_thrust = ac_cfg_mass() * (acc - [0;0;9.80665]);
F_dir = F_thrust / max(norm(F_thrust), 1e-9);
% recover_states does this internally
[F_NED_rec, qd_rec, ~] = df.recover_states(0.5*sum(times), ac_cfg_mass(), 9.80665);
F_dir_rec = F_NED_rec / max(norm(F_NED_rec), 1e-9);
tr = test_helpers('assert_near', tr, F_dir_rec, F_dir, 1e-9, 'T6 thrust direction recovery');

tr = test_helpers('report', tr);
end

function m = ac_cfg_mass()
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'config'));
    cfg = aircraft_config();
    m = cfg.mass;
end
