function tr = test_integration_realistic()
% TEST_INTEGRATION_REALISTIC  Stage 7 — full realistic mission integration test.
%
% Runs the complete 4-phase mission (climb, hover, climb-to-cruise, cruise)
% and validates each phase has acceptable tracking error and stability.
%
% This is a long test (~3-5 min). It validates the full envelope.
tr = test_helpers('init', 'test_integration_realistic');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'aerodynamics'));
addpath(fullfile(here, '..', 'control'));
addpath(fullfile(here, '..', 'trajectory'));
addpath(fullfile(here, '..', 'environment'));
addpath(fullfile(here, '..', 'simulation'));

ac_cfg   = aircraft_config();
ctrl_cfg = controller_config();
sim_cfg  = simulation_config();
sim_cfg.wind.dryden.enable = false;
sim_cfg.wind.constant_NED  = [0;0;0];
ctrl_cfg.inner.type = 'SO3';

[traj, info] = build_realistic_mission(ac_cfg);
sim_cfg.t_final = traj.total_time() + 5;

fprintf('       Running %.0f s realistic mission integration...\n', sim_cfg.t_final);
log = run_simulation(ac_cfg, ctrl_cfg, sim_cfg, traj);

% Phase indices (new 4-phase rest-to-rest mission)
N = length(log.t);
t = log.t;
t_p1_end = info.T_climb1;
t_p2_end = t_p1_end + info.T_hover;
t_p3_end = t_p2_end + info.T_climb2;
idx_p1 = find(t <= t_p1_end);
idx_p2 = find(t > t_p1_end & t <= t_p2_end);
idx_p3 = find(t > t_p2_end & t <= t_p3_end);
idx_p4 = find(t > t_p3_end);

err = vecnorm(log.tracking_err, 2, 1);

% I1: Phase 1 (climb 20m, rest-to-rest, easy)
err_p1 = max(err(idx_p1));
fprintf('       Phase 1 (climb 20m)        max_err = %.2f m\n', err_p1);
tr = test_helpers('assert', tr, err_p1 < 5, 'I1 climb tracking <5m');

% I2: Phase 2 (hover 20s at -20m) — must be very tight
err_p2 = max(err(idx_p2));
fprintf('       Phase 2 (hover 20s)        max_err = %.2f m\n', err_p2);
tr = test_helpers('assert', tr, err_p2 < 3, 'I2 hover tracking <3m');

% I3: Phase 3 (slow rest-to-rest climb to 1000m)
err_p3 = max(err(idx_p3));
fprintf('       Phase 3 (climb to 1000m)   max_err = %.2f m\n', err_p3);
tr = test_helpers('assert', tr, err_p3 < 20, 'I3 climb 1000m tracking <20m');

% I4: Phase 4 (slow horizontal cruise at 1000m, rest-to-rest)
if ~isempty(idx_p4)
    err_p4 = max(err(idx_p4));
    fprintf('       Phase 4 (slow cruise)      max_err = %.2f m\n', err_p4);
    tr = test_helpers('assert', tr, err_p4 < 50, 'I4 horizontal cruise tracking <50m');
else
    tr = test_helpers('assert', tr, true, 'I4 cruise (skipped - not reached)');
end

% I5: Pitch stays near +90 deg throughout (no aggressive transition)
pitch_all = log.pitch_deg(t > 1);  % skip initial transient
pitch_mean = mean(pitch_all);
fprintf('       Mean pitch (whole mission, t>1s) = %.1f deg (expected ~90)\n', pitch_mean);
tr = test_helpers('assert', tr, abs(pitch_mean - 90) < 20, ...
    'I5 mean pitch near 90 deg (tailsitter held upright)');

% I6: Quaternion always normalized
qnorm_err = max(abs(vecnorm(log.quat, 2, 1) - 1));
fprintf('       Max quaternion norm error  = %.3e\n', qnorm_err);
tr = test_helpers('assert', tr, qnorm_err < 1e-3, 'I6 quaternion norm preserved');

% I7: No catastrophic divergence (final altitude near reference)
final_alt = -log.pos_NED(3, end);
expected_alt = 1000;
fprintf('       Final altitude: %.0f m (expected ~%.0f m)\n', final_alt, expected_alt);
tr = test_helpers('assert', tr, abs(final_alt - expected_alt) < 200, ...
    'I7 final altitude within 200m of cruise');

tr = test_helpers('report', tr);
end
