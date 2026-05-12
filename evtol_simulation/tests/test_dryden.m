function tr = test_dryden()
% TEST_DRYDEN  Stage 5b — Dryden turbulence + constant wind.
tr = test_helpers('init', 'test_dryden');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'environment'));

sim = simulation_config();

% W1: disabled returns constant_NED exactly
cfg_disabled = sim.wind;
cfg_disabled.dryden.enable = false;
wm = DrydenWind(cfg_disabled);
v = wm.step(20.0, 100.0, 0.01);
tr = test_helpers('assert_near', tr, v, sim.wind.constant_NED, 1e-15, ...
    'W1 disabled returns constant_NED exactly');

% W2: enabled produces variance > 0
cfg_enabled = sim.wind;
cfg_enabled.dryden.enable = true;
wm2 = DrydenWind(cfg_enabled);
N = 5000;
samples = zeros(3, N);
for k = 1:N
    samples(:, k) = wm2.step(25.0, 50.0, 0.01);
end
% Subtract constant wind
samples = samples - sim.wind.constant_NED;
sigma_emp = std(samples, 0, 2);
fprintf('       W2 empirical sigma: [%.2f %.2f %.2f] m/s\n', sigma_emp);
tr = test_helpers('assert', tr, all(sigma_emp > 0.1), 'W2 turbulence variance > 0');

% W3: reset zeroes states
wm2.reset();
tr = test_helpers('assert_lt', tr, abs(wm2.x_u), 1e-15, 'W3a x_u zeroed');
tr = test_helpers('assert_lt', tr, norm(wm2.x_v), 1e-15, 'W3b x_v zeroed');
tr = test_helpers('assert_lt', tr, norm(wm2.x_w), 1e-15, 'W3c x_w zeroed');

% W4: sample mean ~ zero (turbulent component) — relaxed for finite-sample 3D combined
% With N samples and 3D vector, |mean| ~ sigma * sqrt(3/N_eff) where N_eff is the
% number of independent samples (smaller than N due to filter correlation length).
mean_emp = mean(samples, 2);
mean_err = norm(mean_emp);
% Threshold: 2x max sigma is acceptable for ~50 s of correlated samples
threshold = 2 * max(sigma_emp);
tr = test_helpers('assert_lt', tr, mean_err, threshold, ...
    'W4 turbulent component zero-mean (3D, finite-sample)');

tr = test_helpers('report', tr);
end
