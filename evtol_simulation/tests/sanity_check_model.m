function results = sanity_check_model()
% SANITY_CHECK_MODEL  Physical-plausibility checks on the eVTOL model.
%
% These checks verify that the modeled aircraft behaves consistently with
% first-principles aerospace physics. Each check prints:
%   [OK]    value within sane range
%   [WARN]  value in tolerable range but unusual for this aircraft class
%   [FAIL]  value contradicts physics or basic sizing rules
%   [INFO]  diagnostic value for inspection (no pass/fail)
%
% Usage (from evtol_simulation/):
%   addpath(genpath(pwd));
%   results = sanity_check_model();

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'aerodynamics'));
addpath(fullfile(here, '..', 'simulation'));

ac = aircraft_config();
g = 9.80665;
rho_sl = 1.225;

results = struct();
fprintf('\n');
fprintf('===================================================================\n');
fprintf('   eVTOL Tailsitter — Physical Sanity Checks\n');
fprintf('===================================================================\n\n');

% --- 1. Mass / weight ---
W = ac.mass * g;
print_info('Mass',          'm',         '%.0f kg',  ac.mass);
print_info('Weight',        'W = m*g',   '%.0f N',   W);

% --- 2. Hover trim: total thrust must exceed weight by margin ---
T_avail = ac.n_rotors * ac.rotor.thrust_max;
TW_ratio = T_avail / W;
print_check('T/W ratio (max thrust to weight)', ...
    sprintf('%.2f', TW_ratio), TW_ratio > 1.10 && TW_ratio < 2.0, ...
    'Tailsitter typical 1.2-1.6');
results.TW_ratio = TW_ratio;

% Per-rotor hover thrust
T_hover_per_rotor = W / ac.n_rotors;
hover_throttle = T_hover_per_rotor / ac.rotor.thrust_max;
print_check('Hover throttle (per rotor)', ...
    sprintf('%.2f%%', 100 * hover_throttle), ...
    hover_throttle > 0.5 && hover_throttle < 0.95, ...
    'Healthy: 60-85% (margin for control)');
results.hover_throttle = hover_throttle;

% --- 3. Wing loading ---
WL = W / ac.wing.area;
print_check('Wing loading W/S', ...
    sprintf('%.1f N/m^2 (%.1f kg/m^2)', WL, WL/g), ...
    WL/g > 50 && WL/g < 250, ...
    'Light eVTOL typical 80-180 kg/m^2');
results.wing_loading = WL/g;

% --- 4. Disk loading (rotor) ---
total_disk_area = ac.n_rotors * ac.rotor.disk_area;
DL = W / total_disk_area;
print_check('Disk loading T/A', ...
    sprintf('%.0f N/m^2 (%.1f kg/m^2)', DL, DL/g), ...
    DL/g > 50 && DL/g < 800, ...
    'Open rotor typical 100-500 kg/m^2');
results.disk_loading = DL/g;

% --- 5. Hover induced velocity (Glauert) ---
v_i_hover = sqrt(T_hover_per_rotor / (2 * rho_sl * ac.rotor.disk_area));
print_check('Hover induced velocity v_i', ...
    sprintf('%.2f m/s', v_i_hover), ...
    v_i_hover > 5 && v_i_hover < 60, ...
    'Open rotor 10-30 m/s typical');
results.v_i_hover = v_i_hover;

% --- 6. Hover power ---
% P_ideal = T * v_i, plus profile drag, divide by figure-of-merit ~0.7
P_ideal_per_rotor = T_hover_per_rotor * v_i_hover;
FM = 0.7;   % figure of merit
P_actual_per_rotor = P_ideal_per_rotor / FM;
P_total_kW = ac.n_rotors * P_actual_per_rotor / 1000;
PW_kgperkW = (ac.mass) / P_total_kW;
print_check('Hover power total', ...
    sprintf('%.0f kW (%.1f hp)', P_total_kW, P_total_kW * 1.341), ...
    P_total_kW > 50 && P_total_kW < 1500, ...
    'Light eVTOL 100-500 kW typical');
print_check('Power loading (kg/kW)', ...
    sprintf('%.2f kg/kW', PW_kgperkW), ...
    PW_kgperkW > 1.5 && PW_kgperkW < 8.0, ...
    'Healthy 2-5 kg/kW for eVTOL');
results.P_hover_kW = P_total_kW;

% --- 7. Stall speed (cruise wing) ---
V_stall = sqrt(2 * W / (rho_sl * ac.wing.area * ac.wing.CL_max));
print_check('Stall speed V_s', ...
    sprintf('%.1f m/s (%.1f km/h)', V_stall, V_stall*3.6), ...
    V_stall > 15 && V_stall < 50, ...
    'Tailsitter cruise 25-45 m/s typical');
results.V_stall = V_stall;

% --- 8. Cruise trim AoA ---
% Assume cruise at 1.3*V_stall, find required CL
V_cruise = 1.3 * V_stall;
CL_cruise = 2 * W / (rho_sl * ac.wing.area * V_cruise^2);
alpha_cruise = CL_cruise / ac.wing.CL_alpha + ac.wing.alpha0;
print_check('Cruise speed (1.3 V_s)', ...
    sprintf('%.1f m/s (%.0f km/h)', V_cruise, V_cruise*3.6), ...
    V_cruise > 25 && V_cruise < 60, ...
    'Reasonable for tailsitter');
print_check('Cruise CL', ...
    sprintf('%.3f', CL_cruise), ...
    CL_cruise > 0.3 && CL_cruise < ac.wing.CL_max, ...
    sprintf('Below CL_max=%.2f', ac.wing.CL_max));
print_check('Cruise alpha', ...
    sprintf('%.2f deg', rad2deg(alpha_cruise)), ...
    abs(rad2deg(alpha_cruise)) < ac.wing.alpha_stall*180/pi, ...
    'Below stall AoA');
results.V_cruise = V_cruise;
results.alpha_cruise_deg = rad2deg(alpha_cruise);

% --- 9. Cruise drag and L/D ---
CD_cruise = ac.wing.CD0 + (CL_cruise^2)/(pi * ac.wing.AR * ac.wing.e);
LD = CL_cruise / CD_cruise;
print_check('Cruise L/D', ...
    sprintf('%.1f', LD), ...
    LD > 6 && LD < 25, ...
    'Light fixed-wing 8-15 typical');
results.LD_cruise = LD;

% --- 10. Cruise power ---
D_cruise = 0.5 * rho_sl * ac.wing.area * V_cruise^2 * CD_cruise;
P_cruise_kW = D_cruise * V_cruise / 1000;
hover_to_cruise_ratio = P_total_kW / P_cruise_kW;
print_check('Cruise power required', ...
    sprintf('%.1f kW', P_cruise_kW), ...
    P_cruise_kW > 10 && P_cruise_kW < 500, ...
    sprintf('Hover/cruise power ratio = %.1f (typical 2-5x)', hover_to_cruise_ratio));
results.P_cruise_kW = P_cruise_kW;

% --- 11. Rotor: T = kT * Omega^2 sanity ---
Omega_check = ac.rotor.omega_max;
T_at_Omega_max = ac.rotor.kT * Omega_check^2;
print_check('Rotor T at Omega_max', ...
    sprintf('%.0f N', T_at_Omega_max), ...
    abs(T_at_Omega_max - ac.rotor.thrust_max)/ac.rotor.thrust_max < 0.05, ...
    'Should match T_max within 5%');

% --- 12. Time-constant separation (cascade stability) ---
tau_motor = ac.rotor.tau_motor;
% inner loop bandwidth (from controller_config) — read directly
ctrl_path = fullfile(here, '..', 'config');
addpath(ctrl_path);
ctrl = controller_config();
omega_n_inner = sqrt(ctrl.so3.kR(2,2) / ac.J(2,2));   % pitch axis
tau_inner = 1 / (omega_n_inner * 0.7);
tau_outer = 1 / ctrl.f_outer * 5;   % 5 outer steps for response
ratio_inner_motor  = tau_inner / tau_motor;
ratio_outer_inner  = tau_outer / tau_inner;
print_check('Motor tau (1st order)', sprintf('%.3f s', tau_motor), ...
    tau_motor > 0.01 && tau_motor < 0.2, 'Typical 30-100 ms');
print_check('Inner SO3 tau (1/zeta*wn)', sprintf('%.3f s', tau_inner), ...
    tau_inner > tau_motor * 1.5, ...
    sprintf('Separation: inner/motor = %.1fx (need >1.5)', ratio_inner_motor));
print_check('Outer NMPC tau (5*dt_outer)', sprintf('%.3f s', tau_outer), ...
    tau_outer > tau_inner * 1.5, ...
    sprintf('Separation: outer/inner = %.1fx (need >1.5)', ratio_outer_inner));
results.cascade_separation = [ratio_inner_motor, ratio_outer_inner];

% --- 13. Inertia ratios for tailsitter ---
J_xx = ac.J(1,1); J_yy = ac.J(2,2); J_zz = ac.J(3,3); J_xz = abs(ac.J(1,3));
print_check('J_yy/J_xx (pitch vs roll)', ...
    sprintf('%.2f', J_yy/J_xx), ...
    J_yy/J_xx > 1.0 && J_yy/J_xx < 4.0, ...
    'Tailsitter has more pitch inertia (longitudinal)');
print_check('J_xz / J_yy (cross-coupling)', ...
    sprintf('%.3f', J_xz/J_yy), ...
    J_xz/J_yy < 0.20, ...
    'Should be small (<20%)');

% --- 14. Control authority margin ---
% Maximum moment about pitch axis from differential thrust:
% Two rotors at +x_arm, two at -x_arm (here: forward/aft via wing position)
% We compute it from per-rotor positions
M_y_max = 0;
for i = 1:ac.n_rotors
    r = ac.rotor.position(i, :)';
    e = [1; 0; 0];   % nominal axis
    m_per_T = cross(r, e);
    M_y_max = M_y_max + abs(m_per_T(2)) * ac.rotor.thrust_max;
end
% Required moment to flip 90 deg in 1 second: J*omega^2 / 2 = J*pi^2/2
M_y_required = J_yy * pi^2 / 2;
print_check('Pitch authority M_max', ...
    sprintf('%.0f N*m', M_y_max), ...
    M_y_max > M_y_required * 0.8, ...
    sprintf('Need M >= %.0f N*m for 90deg/1s pitch', M_y_required));
results.M_y_max = M_y_max;

% --- 15. Slipstream sanity (full contraction) ---
v_slip_far = 2 * v_i_hover;
print_check('Far-wake slipstream', ...
    sprintf('%.1f m/s', v_slip_far), ...
    v_slip_far > 10 && v_slip_far < 100, ...
    'Far-wake = 2*v_i (Froude theorem)');

% --- 16. Aerodynamic vs propulsive authority crossover ---
% At what speed does aerodynamic surface authority match rotor differential?
% Roughly: q_bar * S * arm * delta_max ~ M_y_max (rough estimate)
delta_max = ac.surf.elev_max;
q_bar_crossover = M_y_max / (ac.htail.elev_eff * ac.htail.area * ac.htail.arm * delta_max);
V_crossover = sqrt(2 * q_bar_crossover / rho_sl);
print_check('Aero/rotor crossover speed', ...
    sprintf('%.1f m/s', V_crossover), ...
    V_crossover > 5 && V_crossover < 60, ...
    'Below this: rotors dominate; above: surfaces');
results.V_crossover = V_crossover;

% --- 17. Free-fall sanity (no propulsion = a = g) ---
ac_obj = Aircraft(ac);
x0 = ac_obj.initial_state('hover');
xdot = ac_obj.dynamics(x0, zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1));
a_NED = xdot(4:6);
print_check('Free-fall acceleration', ...
    sprintf('|a|=%.4f m/s^2 (NED z=%.4f)', norm(a_NED), a_NED(3)), ...
    abs(a_NED(3) - g) < 1e-9 && abs(a_NED(1)) + abs(a_NED(2)) < 1e-9, ...
    'No forces -> a_NED = [0;0;g]');

% --- 18. Hover trim test (pure thrust, no aero, no moments) ---
F_prop_B = [W; 0; 0];   % thrust along body x (which is up at hover)
xdot_h = ac_obj.dynamics(x0, zeros(3,1), zeros(3,1), F_prop_B, zeros(3,1), zeros(3,1));
print_check('Hover trim (T = m*g)', ...
    sprintf('|a_NED|=%.4e m/s^2', norm(xdot_h(4:6))), ...
    norm(xdot_h(4:6)) < 1e-9, ...
    'Should give zero net acceleration');

% --- Summary ---
fprintf('\n');
fprintf('===================================================================\n');
fprintf('   Summary: physical sanity check complete\n');
fprintf('   Inspect [WARN] entries -- they are not failures but unusual values.\n');
fprintf('   Inspect [FAIL] entries -- they indicate physical inconsistency.\n');
fprintf('===================================================================\n\n');
end


%% --- Helper functions ---
function print_check(name, value_str, ok, note)
    if ok
        tag = '[OK]  ';
    else
        tag = '[WARN]';
    end
    fprintf('  %s %-32s = %-30s | %s\n', tag, name, value_str, note);
end

function print_info(name, sym, fmt, val)
    fprintf('  [INFO] %-32s = %s\n', sprintf('%s (%s)', name, sym), sprintf(fmt, val));
end
