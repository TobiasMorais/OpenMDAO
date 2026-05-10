function tr = test_force_to_attitude()
% TEST_FORCE_TO_ATTITUDE  Stage 3e — F_des -> (q_d, T_total, omega_d) mapping.
% This is the PRIMARY suspect for the cascaded divergence observed in main run.
tr = test_helpers('init', 'test_force_to_attitude');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'control'));

ac_cfg = aircraft_config();
g = 9.80665;
m = ac_cfg.mass;

% F1: hover specific thrust = -g_NED -> qd should give x_B = up (NED -Z)
% NMPC outputs f_specific = -g_NED for hover (equiv to F_NED = -m g_NED)
F_des_NED = [0; 0; -g];   % specific thrust pointing up (NED -Z)
[qd, T_total, Wd, Wd_dot] = force_to_attitude(F_des_NED, 0.0, 0.0, []);
R_d = quat_utils('toR', qd);
x_B_world = R_d * [1;0;0];
fprintf('       F1 hover: x_B in world = [%.3f %.3f %.3f]\n', x_B_world);
% expected x_B = [0; 0; -1] (up)
tr = test_helpers('assert_near', tr, x_B_world, [0;0;-1], 1e-6, ...
    'F1 hover: x_B points up in NED');

% F2: forward acceleration -> qd tilts forward
a_fwd = 5.0;
F_des_NED = [m*a_fwd; 0; -m*g] / m;   % specific thrust
F_des_NED = [a_fwd; 0; -g];
[qd, ~, ~, ~] = force_to_attitude(F_des_NED, 0.0, 0.0, []);
R_d = quat_utils('toR', qd);
x_B_world = R_d * [1;0;0];
% x_B should be parallel to F_des direction
F_dir = F_des_NED / norm(F_des_NED);
dot_prod = x_B_world' * F_dir;
fprintf('       F2 fwd accel: x_B.F_dir = %.6f (expect 1)\n', dot_prod);
tr = test_helpers('assert', tr, abs(dot_prod - 1) < 1e-3, 'F2 x_B aligned with F_des');

% F3: continuity (Lipschitz) — small dF -> small dqd
F_a = [0.1; 0; -g];
F_b = [0.11; 0; -g];
[qd_a, ~,~,~] = force_to_attitude(F_a, 0.0, 0.0, []);
[qd_b, ~,~,~] = force_to_attitude(F_b, 0.0, 0.0, []);
qe = quat_utils('errMul', qd_a, qd_b);
ang_err = 2 * acos(max(-1, min(1, qe(1))));
fprintf('       F3 continuity: small dF=0.01, ang_err=%.3e rad\n', ang_err);
tr = test_helpers('assert_lt', tr, ang_err, 0.01, 'F3 continuity (Lipschitz)');

% F4: det(R(qd)) = +1
[qd, ~,~,~] = force_to_attitude([0; 0.5; -g], 0.0, 0.0, []);
R_d = quat_utils('toR', qd);
tr = test_helpers('assert', tr, abs(det(R_d) - 1) < 1e-6, 'F4 det(R(qd)) = +1');

% F5: x_B parallel to F_des in general direction
F_test = [3; 1; -8];
[qd, ~,~,~] = force_to_attitude(F_test, 0.0, 0.0, []);
R_d = quat_utils('toR', qd);
x_B = R_d * [1;0;0];
F_dir = F_test / norm(F_test);
fprintf('       F5 generic: x_B.F_dir = %.6f\n', x_B'*F_dir);
tr = test_helpers('assert', tr, x_B'*F_dir > 0.999, 'F5 x_B || F_des (generic)');

% F6: heading rotation -> y_B varies in plane perp to x_B
F_const = [0; 0; -g];
[qd_psi0, ~,~,~] = force_to_attitude(F_const, 0.0, 0.0, []);
[qd_psi1, ~,~,~] = force_to_attitude(F_const, deg2rad(45), 0.0, []);
y_B_0 = quat_utils('toR', qd_psi0) * [0;1;0];
y_B_1 = quat_utils('toR', qd_psi1) * [0;1;0];
fprintf('       F6 heading: y_B(0)=[%.2f %.2f %.2f] y_B(45)=[%.2f %.2f %.2f]\n', ...
    y_B_0, y_B_1);
% Both must be perpendicular to x_B = [0;0;-1] (z-component must be ~0 for hover)
tr = test_helpers('assert_lt', tr, abs(y_B_0(3))+abs(y_B_1(3)), 0.05, ...
    'F6a y_B perp to x_B');
% But y_B should have rotated
diff_yB = norm(y_B_0 - y_B_1);
tr = test_helpers('assert', tr, diff_yB > 0.1, 'F6b heading change moves y_B');

% F7: no NaN at hover degenerate (force exactly along -Z)
F_deg = [0; 0; -m*g];   % large but degenerate hover
[qd, T, Wd, Wd_dot] = force_to_attitude(F_deg, 0.0, 0.0, []);
tr = test_helpers('assert', tr, all(isfinite(qd)) && all(isfinite(Wd)), ...
    'F7 no NaN at degenerate hover');

tr = test_helpers('report', tr);
end
