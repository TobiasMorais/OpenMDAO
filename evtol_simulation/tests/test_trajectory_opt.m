function tr = test_trajectory_opt()
% TEST_TRAJECTORY_OPT  Stage 4b — Trajectory optimizer (Flat baseline + GA).
tr = test_helpers('init', 'test_trajectory_opt');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'trajectory'));

ac_cfg = aircraft_config();
wp = [ 0, 0,   0, 0;
       0, 0, -50, 0;
     100, 0, -50, 0;
     150, 0, -10, 0];
opt = TrajectoryOptimizer(wp, ac_cfg);

% O1: baseline total time consistent
traj = opt.optimize('method', 'flat');
T_total = traj.total_time();
tr = test_helpers('assert', tr, T_total > 0 && T_total < 100, 'O1 baseline total time sane');

% O2: GA cost <= baseline cost (only if Global Optimization Toolbox available)
has_ga = exist('ga', 'file') == 2;
if has_ga
    fprintf('       Running GA (this may take ~10-30 s)...\n');
    J_base = opt.energy_cost(ones(size(wp,1)-1, 1), opt.baseline.time_alloc);
    traj_ga = opt.optimize('method', 'ga', 'ga_pop', 16, 'ga_gens', 6);
    s_opt = traj_ga.time_alloc ./ opt.baseline.time_alloc;
    J_ga = opt.energy_cost(s_opt, opt.baseline.time_alloc);
    fprintf('       J_baseline=%.3e, J_ga=%.3e\n', J_base, J_ga);
    tr = test_helpers('assert', tr, J_ga <= J_base * 1.05, 'O2 GA does not increase cost');
else
    fprintf('       (Global Optimization Toolbox not present; skipping O2)\n');
    tr = test_helpers('assert', tr, true, 'O2 GA skipped (toolbox missing)');
end

% O3: no thrust saturation in baseline
N_samples = 50;
ts = linspace(0, T_total, N_samples);
saturated = false;
for k = 1:N_samples
    [~, ~, acc, ~, ~, ~, ~] = traj.eval(ts(k));
    F_NED = ac_cfg.mass * (acc - [0;0;9.80665]);
    T_per_rotor = norm(F_NED) / ac_cfg.n_rotors;
    if T_per_rotor > ac_cfg.rotor.thrust_max
        saturated = true;
        break;
    end
end
tr = test_helpers('assert', tr, ~saturated, 'O3 baseline trajectory feasible (no saturation)');

tr = test_helpers('report', tr);
end
