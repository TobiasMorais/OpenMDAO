function tr = test_so3_controller()
% TEST_SO3_CONTROLLER  Stage 3a — geometric SO(3) attitude controller.
% This suite is EXPECTED to surface the gain-bandwidth issue (C3) for
% the 1500 kg vehicle: gains tuned for quadrotors are too small here.
tr = test_helpers('init', 'test_so3_controller');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'control'));

ac_cfg   = aircraft_config();
ctrl_cfg = controller_config();
ctrl     = AttitudeControllerSO3(ctrl_cfg);

% C1: equilibrium R = Rd, omega = 0
R  = eye(3);
Rd = eye(3);
omega = [0;0;0];
M = ctrl.compute(R, omega, Rd, [0;0;0], [0;0;0], ac_cfg.J, 0);
tr = test_helpers('assert_lt', tr, norm(M), 1e-9, 'C1 equilibrium => M = 0');

% C2: small attitude error in y (pitch) -> M_y should counter (negative)
ctrl.reset();
Rd = so3_utils('exp', [0; 0.01; 0]);
M2 = ctrl.compute(eye(3), [0;0;0], Rd, [0;0;0], [0;0;0], ac_cfg.J, 0);
% e_R = 0.5 * vee(Rd^T R - R^T Rd). For small theta_y: Rd ~ I + hat(0.01 ey), so
% Rd^T R - R^T Rd = -hat(0.01 ey) - hat(0.01 ey) = -2 hat(0.01 ey)
% vee(-2 hat(0.01 ey)) = [0; -0.02; 0], 0.5 * = [0; -0.01; 0]
% M = -kR * e_R = -diag(...) * [0; -0.01; 0] = [0; +0.01*kR_yy; 0] -> POSITIVE
tr = test_helpers('assert', tr, M2(2) > 0, 'C2 sign of restoring moment correct');

% C3: BANDWIDTH CHECK — omega_n = sqrt(kR/J) >= 5 rad/s for tailsitter
J_yy = ac_cfg.J(2,2);
J_xx = ac_cfg.J(1,1);
J_zz = ac_cfg.J(3,3);
kR_xx = ctrl_cfg.so3.kR(1,1);
kR_yy = ctrl_cfg.so3.kR(2,2);
kR_zz = ctrl_cfg.so3.kR(3,3);
omega_n_x = sqrt(kR_xx / J_xx);
omega_n_y = sqrt(kR_yy / J_yy);
omega_n_z = sqrt(kR_zz / J_zz);
min_bw = min([omega_n_x, omega_n_y, omega_n_z]);
fprintf('       Bandwidth (rad/s): x=%.2f y=%.2f z=%.2f\n', omega_n_x, omega_n_y, omega_n_z);
tr = test_helpers('assert', tr, min_bw >= 5.0, 'C3 attitude bandwidth >= 5 rad/s', ...
    sprintf('min omega_n=%.2f rad/s (need >=5 for tailsitter)', min_bw));

% C4: damping ratio in [0.5, 1.5]
zeta_x = ctrl_cfg.so3.kOm(1,1) / (2 * sqrt(kR_xx * J_xx));
zeta_y = ctrl_cfg.so3.kOm(2,2) / (2 * sqrt(kR_yy * J_yy));
zeta_z = ctrl_cfg.so3.kOm(3,3) / (2 * sqrt(kR_zz * J_zz));
fprintf('       Damping zeta: x=%.2f y=%.2f z=%.2f\n', zeta_x, zeta_y, zeta_z);
tr = test_helpers('assert', tr, zeta_x>0.3 && zeta_y>0.3 && zeta_z>0.3 && zeta_x<2 && zeta_y<2 && zeta_z<2, ...
    'C4 damping ratio in reasonable band [0.3, 2.0]');

% C5: integral anti-windup
ctrl.reset();
Rd_persist = so3_utils('exp', [0; 0.5; 0]);
for t = 0:0.01:5
    ctrl.compute(eye(3), [0;0;0], Rd_persist, [0;0;0], [0;0;0], ac_cfg.J, t);
end
tr = test_helpers('assert_lt', tr, max(abs(ctrl.ei)), ctrl_cfg.so3.I_max + 1e-6, ...
    'C5 integral saturated to I_max');

% C6: simulated step recovery
ctrl.reset();
J = ac_cfg.J;
R = eye(3);
omega = [0;0;0];
Rd = so3_utils('exp', [0; deg2rad(30); 0]);
dt = 0.002;
T_settle = NaN;
for k = 1:5000
    t = k*dt;
    M = ctrl.compute(R, omega, Rd, [0;0;0], [0;0;0], J, t);
    omega = omega + dt * (J \ M);
    R = R * so3_utils('exp', omega * dt);
    e_R = so3_utils('errMat', R, Rd);
    if norm(e_R) < 0.05 && isnan(T_settle)
        T_settle = t;
        break;
    end
end
fprintf('       Settling time (5%% on 30deg step) = %.2f s\n', T_settle);
tr = test_helpers('assert', tr, ~isnan(T_settle) && T_settle < 5.0, ...
    'C6 30deg step settles in <5s', sprintf('T_settle=%.2fs', T_settle));

tr = test_helpers('report', tr);
end
