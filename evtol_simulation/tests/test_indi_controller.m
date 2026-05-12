function tr = test_indi_controller()
% TEST_INDI_CONTROLLER  Stage 3b — INDI attitude controller.
tr = test_helpers('init', 'test_indi_controller');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'control'));

ac_cfg   = aircraft_config();
ctrl_cfg = controller_config();
ac       = Aircraft(ac_cfg);
alloc    = ControlAllocator(ctrl_cfg, ac_cfg, ac.rotor_axes);
G_eff    = alloc.attitude_effectiveness();
indi     = AttitudeControllerINDI(ctrl_cfg, G_eff, ac_cfg.n_rotors);

% N1: reset zeros all filter states
indi.reset();
tr = test_helpers('assert_lt', tr, norm(indi.omega_dot_filt), 1e-15, 'N1 reset omega_dot_filt');
tr = test_helpers('assert_lt', tr, norm(indi.u_filt), 1e-15, 'N1 reset u_filt');

% N2: low-pass converges on step input
indi.reset();
omega_traj = zeros(3, 200);
omega_traj(1, 50:end) = 1.0;   % step in p
dt = 0.01;
omega_dot_filt_history = zeros(1, 200);
for k = 1:200
    t = k*dt;
    if k>1
        % manually drive prev_omega for FD
        indi.compute([1;0;0;0], omega_traj(:,k), [1;0;0;0], [0;0;0], zeros(4,1), t);
    end
    omega_dot_filt_history(k) = indi.omega_dot_filt(1);
end
% After ~5 time-constants (1/lp_omega = 0.02s), filtered value should be small (input went constant)
tr = test_helpers('assert', tr, max(omega_dot_filt_history(end-20:end)) < ...
    max(omega_dot_filt_history) * 0.5, 'N2 LP converges after step');

% N3: at equilibrium q=qd, omega=0 -> Du small
indi.reset();
[Du, ~] = indi.compute([1;0;0;0], [0;0;0], [1;0;0;0], [0;0;0], zeros(4,1), 0);
tr = test_helpers('assert_lt', tr, norm(Du), 1e-6, 'N3 equilibrium => Du ~ 0');

% N4: G^+ well-conditioned (singular value ratio)
G = G_eff;
sv = svd(G);
cond_ratio = sv(1) / sv(end);
fprintf('       G condition: smax=%.3e, smin=%.3e, ratio=%.2e\n', sv(1), sv(end), cond_ratio);
tr = test_helpers('assert', tr, cond_ratio < 1e6, 'N4 G well-conditioned', ...
    sprintf('cond=%.2e', cond_ratio));

% N5: response to qd != q produces nonzero Du
indi.reset();
qd = quat_utils('fromEuler', 0, 0.1, 0);  % small pitch
[Du, ~] = indi.compute([1;0;0;0], [0;0;0], qd, [0;0;0], zeros(4,1), 0);
tr = test_helpers('assert', tr, norm(Du) > 1e-6, 'N5 qd != q => Du nonzero');

tr = test_helpers('report', tr);
end
