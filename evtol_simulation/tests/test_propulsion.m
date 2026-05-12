function tr = test_propulsion()
% TEST_PROPULSION  Stage 2b — rotor model: BEMT + ESC dynamics + induced velocity.
tr = test_helpers('init', 'test_propulsion');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'aerodynamics'));

ac_cfg = aircraft_config();
ac     = Aircraft(ac_cfg);
prop   = Propulsion(ac_cfg, ac.rotor_axes);

% P1: zero rotor speed => zero thrust
[F0, M0, vi0, ~] = prop.compute([0;0;0], 1.225);
tr = test_helpers('assert_lt', tr, norm(F0)+norm(M0), 1e-9, 'P1 Omega=0 => T=0');

% P2: T = kT * Omega^2 (static)
prop.Omega_actual = 0.5 * ac_cfg.rotor.omega_max * ones(4,1);
[F2, ~, ~, ~] = prop.compute([0;0;0], 1.225);
expected_T_per = ac_cfg.rotor.kT * (0.5*ac_cfg.rotor.omega_max)^2;
% Net thrust along body x is sum of axial components (cants are small)
T_total_expected = 4 * expected_T_per;
% F_prop_B is along rotor_axes (mostly +X_B); take dot product
T_along_x = F2(1);
tr = test_helpers('assert', tr, abs(T_along_x - T_total_expected)/T_total_expected < 0.02, ...
    'P2 static T = kT*Omega^2', sprintf('actual=%.0f expected=%.0f', T_along_x, T_total_expected));

% P3: ESC first-order convergence — 5*tau should give ~99% of step
prop = Propulsion(ac_cfg, ac.rotor_axes);
T_cmd = 2000 * ones(4,1);
dt = 0.001;
N = round(5 * ac_cfg.rotor.tau_motor / dt);
for k = 1:N
    prop.step_actuator(T_cmd, dt);
end
T_actual = ac_cfg.rotor.kT * prop.Omega_actual.^2;
tr = test_helpers('assert', tr, all(T_actual > 0.98 * T_cmd), 'P3 ESC reaches 98% in 5*tau', ...
    sprintf('min T_act/T_cmd = %.3f', min(T_actual./T_cmd)));

% P4: saturation
prop = Propulsion(ac_cfg, ac.rotor_axes);
T_huge = 1e5 * ones(4,1);
for k = 1:100
    prop.step_actuator(T_huge, dt);
end
T_sat = ac_cfg.rotor.kT * prop.Omega_actual.^2;
tr = test_helpers('assert', tr, all(T_sat <= ac_cfg.rotor.thrust_max + 1), 'P4 saturation at T_max', ...
    sprintf('max T = %.0f (limit %.0f)', max(T_sat), ac_cfg.rotor.thrust_max));

% P5: induced velocity hover
prop = Propulsion(ac_cfg, ac.rotor_axes);
T_hover = ac_cfg.mass * ac.gravity / 4;     % per rotor
prop.Omega_actual = sqrt(T_hover / ac_cfg.rotor.kT) * ones(4,1);
[~, ~, vi, ~] = prop.compute([0;0;0], 1.225);
v_i_expected = sqrt(T_hover / (2 * 1.225 * ac_cfg.rotor.disk_area));
tr = test_helpers('assert', tr, abs(vi(1) - v_i_expected)/v_i_expected < 0.02, ...
    'P5 hover v_i = sqrt(T/(2 rho A))', sprintf('vi=%.2f vexp=%.2f', vi(1), v_i_expected));

% P6: forward flight reduces v_i
[~, ~, vi_fwd, ~] = prop.compute([30; 0; 0], 1.225);
tr = test_helpers('assert', tr, vi_fwd(1) < vi(1), 'P6 v_i decreases with V_x', ...
    sprintf('hover %.2f, fwd %.2f', vi(1), vi_fwd(1)));

% P7: reaction torque sign
prop = Propulsion(ac_cfg, ac.rotor_axes);
prop.Omega_actual = 100 * ones(4,1);
[~, M, ~, ~] = prop.compute([0;0;0], 1.225);
% with paired CCW/CW reaction torques cancel; check magnitude
tr = test_helpers('assert_lt', tr, abs(M(1)) + abs(M(2)) + abs(M(3)), 100, ...
    'P7 paired reactions partially cancel');

% P8: hover trim equality
prop = Propulsion(ac_cfg, ac.rotor_axes);
T_per = ac_cfg.mass * ac.gravity / 4;
prop.Omega_actual = sqrt(T_per / ac_cfg.rotor.kT) * ones(4,1);
[Fp, ~, ~, ~] = prop.compute([0;0;0], 1.225);
tr = test_helpers('assert', tr, abs(Fp(1) - ac_cfg.mass*ac.gravity)/abs(ac_cfg.mass*ac.gravity) < 0.05, ...
    'P8 hover equality: 4 rotors balance gravity', ...
    sprintf('F_x=%.0f required=%.0f', Fp(1), ac_cfg.mass*ac.gravity));

tr = test_helpers('report', tr);
end
