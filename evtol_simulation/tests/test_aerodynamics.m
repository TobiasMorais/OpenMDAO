function tr = test_aerodynamics()
% TEST_AERODYNAMICS  Stage 2a — Viterna-Corrigan + strip theory.
tr = test_helpers('init', 'test_aerodynamics');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'aerodynamics'));

ac_cfg = aircraft_config();
aero   = Aerodynamics(ac_cfg);
defs   = aero.build_surface_strips();

% A1: zero airspeed => zero force
[F0, M0] = aero.compute_forces_moments([0;0;0], [0;0;0], defs, zeros(3, numel(defs)), 1.225);
tr = test_helpers('assert_lt', tr, norm(F0)+norm(M0), 1e-6, 'A1 V=0 => F=0,M=0');

% A2: NACA-0015 at alpha=0 => CL~0
[CL0, CD0] = aero.airfoil_360(0);
tr = test_helpers('assert_lt', tr, abs(CL0), 0.05, 'A2 alpha=0 => CL ~ 0');
tr = test_helpers('assert_lt', tr, abs(CD0 - ac_cfg.wing.CD0), 0.005, 'A2b alpha=0 => CD ~ CD0');

% A3: linear regime CL = CL_alpha * alpha
alphas = deg2rad([-5, -2, 0, 2, 5, 8]);
errs = zeros(size(alphas));
for k = 1:numel(alphas)
    [CLk, ~] = aero.airfoil_360(alphas(k));
    expected = ac_cfg.wing.CL_alpha * alphas(k);
    errs(k) = abs(CLk - expected);
end
tr = test_helpers('assert_lt', tr, max(errs), 0.05, 'A3 linear CL regime');

% A4: post-stall CL drop
[CL_15, ~] = aero.airfoil_360(deg2rad(15));
[CL_20, ~] = aero.airfoil_360(deg2rad(20));
tr = test_helpers('assert', tr, abs(CL_20) < abs(CL_15), 'A4 post-stall CL decrease', ...
    sprintf('CL(15)=%.3f CL(20)=%.3f', CL_15, CL_20));

% A5: flat plate at 90deg
[CL_90, CD_90] = aero.airfoil_360(deg2rad(90));
tr = test_helpers('assert', tr, abs(CD_90) > 0.8 && abs(CD_90) < 1.5, ...
    'A5 CD(90deg) ~ flat plate', sprintf('CD=%.2f', CD_90));

% A6: symmetric airfoil antisymmetry
[CLp, ~] = aero.airfoil_360(deg2rad(8));
[CLn, ~] = aero.airfoil_360(deg2rad(-8));
tr = test_helpers('assert_lt', tr, abs(CLp + CLn), 0.05, 'A6 symmetric: CL(-a) = -CL(a)');

% A7: dynamic pressure scaling
[F1, ~] = aero.compute_forces_moments([20;0;0], [0;0;0], defs, zeros(3, numel(defs)), 1.225);
[F2, ~] = aero.compute_forces_moments([40;0;0], [0;0;0], defs, zeros(3, numel(defs)), 1.225);
ratio = norm(F2) / max(norm(F1), 1e-9);
tr = test_helpers('assert', tr, abs(ratio - 4) < 0.5, ...
    'A7 F scales with V^2', sprintf('ratio=%.2f (expected 4)', ratio));

% A8: slipstream raises V_local on bathed strip
V_inf_zero = [0;0;0];
V_slip = zeros(3, numel(defs));
V_slip(:, 1) = [25; 0; 0];   % wing_R bathed by rotor wash
[F_no, ~] = aero.compute_forces_moments(V_inf_zero, [0;0;0], defs, zeros(3,numel(defs)), 1.225);
[F_yes, ~] = aero.compute_forces_moments(V_inf_zero, [0;0;0], defs, V_slip, 1.225);
tr = test_helpers('assert', tr, norm(F_yes) > norm(F_no), 'A8 slipstream produces force', ...
    sprintf('F_no=%.2f F_yes=%.2f', norm(F_no), norm(F_yes)));

tr = test_helpers('report', tr);
end
