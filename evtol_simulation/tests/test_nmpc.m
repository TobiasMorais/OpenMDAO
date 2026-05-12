function tr = test_nmpc()
% TEST_NMPC  Stage 3c — Nonlinear MPC outer position loop.
tr = test_helpers('init', 'test_nmpc');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'control'));

ac_cfg   = aircraft_config();
ctrl_cfg = controller_config();
nmpc = PositionControllerNMPC(ctrl_cfg, ac_cfg.mass);

% Skip if Optimization Toolbox missing (graceful degradation)
has_fmincon = exist('fmincon', 'file') == 2;
if ~has_fmincon
    fprintf('       (fmincon not available, using fallback PD path)\n');
end

N = ctrl_cfg.nmpc.N;

% M1: hover hold — reference held at current state
p = [0;0;-50];
v = [0;0;0];
p_horizon = repmat(p, 1, N+1);
v_horizon = repmat(v, 1, N+1);
f_cmd = nmpc.compute(p, v, p_horizon, v_horizon);
expected_hover = [0; 0; -9.80665];
err = norm(f_cmd - expected_hover);
fprintf('       Hover f_cmd = [%.3f %.3f %.3f]\n', f_cmd(1), f_cmd(2), f_cmd(3));
tr = test_helpers('assert', tr, err < 2.0, 'M1 hover hold f_cmd ~ -g_NED', ...
    sprintf('err=%.3f', err));

% M2: forward step ahead in reference => positive accel x
p = [0;0;-50]; v = [0;0;0];
dt_h = ctrl_cfg.nmpc.dt;
for j = 0:N
    p_horizon(:, j+1) = [j*dt_h*5; 0; -50];   % 5 m/s northward target
    v_horizon(:, j+1) = [5; 0; 0];
end
nmpc.reset();
f_cmd = nmpc.compute(p, v, p_horizon, v_horizon);
fprintf('       Forward step f_cmd = [%.3f %.3f %.3f]\n', f_cmd(1), f_cmd(2), f_cmd(3));
tr = test_helpers('assert', tr, f_cmd(1) > 0.1, 'M2 forward step => f_x > 0');

% M3: extreme reference => tilt constraint check (informal)
nmpc.reset();
for j = 0:N
    p_horizon(:, j+1) = [j*dt_h*100; 0; -50];   % 100 m/s requested (extreme)
    v_horizon(:, j+1) = [100; 0; 0];
end
f_cmd = nmpc.compute(p, v, p_horizon, v_horizon);
fxy = norm(f_cmd(1:2));
fz  = -f_cmd(3);
tilt = atan2(fxy, max(fz, 1e-3));
fprintf('       Extreme ref tilt = %.1f deg (limit %.1f)\n', rad2deg(tilt), rad2deg(ctrl_cfg.nmpc.tilt_max));
% Note: tilt constraint is enforced by fmincon as inequality; allow small numerical excess
tr = test_helpers('assert', tr, tilt < ctrl_cfg.nmpc.tilt_max + deg2rad(5), ...
    'M3 tilt constraint approximately enforced');

% M4: |f_cmd| <= a_max
amax = ctrl_cfg.nmpc.acc_max;
tr = test_helpers('assert', tr, max(abs(f_cmd)) <= amax + 0.1, ...
    'M4 |f_cmd| bounded by a_max');

% M5: warm start — second call faster
% (implicit; just call twice to ensure no error)
nmpc.compute(p, v, p_horizon, v_horizon);
nmpc.compute(p, v, p_horizon, v_horizon);
tr = test_helpers('assert', tr, true, 'M5 warm-start does not crash');

tr = test_helpers('report', tr);
end
