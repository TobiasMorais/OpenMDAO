function tr = test_atmosphere()
% TEST_ATMOSPHERE  Stage 5c — ISA standard atmosphere density.
tr = test_helpers('init', 'test_atmosphere');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'environment'));
addpath(fullfile(here, '..', 'config'));

sim = simulation_config();

% Z1: rho(0) = 1.225 kg/m^3
rho_sl = atmosphere_isa(0, sim.env);
tr = test_helpers('assert_lt', tr, abs(rho_sl - 1.225), 1e-3, 'Z1 rho(0) = 1.225 kg/m^3');

% Z2: rho(11000) reasonable
rho_top = atmosphere_isa(11000, sim.env);
fprintf('       rho(11 km) = %.4f kg/m^3 (expected ~0.36)\n', rho_top);
tr = test_helpers('assert', tr, rho_top > 0.3 && rho_top < 0.45, ...
    'Z2 rho(11km) in tropopause range');

% Z3: monotone decrease
heights = [0, 500, 1500, 3000, 6000, 10000];
rhos = arrayfun(@(h) atmosphere_isa(h, sim.env), heights);
diffs = diff(rhos);
tr = test_helpers('assert', tr, all(diffs < 0), 'Z3 rho monotone decreasing');

tr = test_helpers('report', tr);
end
