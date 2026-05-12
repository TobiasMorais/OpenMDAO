function [traj, info] = build_realistic_mission(ac_cfg, opts)
% BUILD_REALISTIC_MISSION  Conservative tailsitter mission (rest-to-rest only).
%
% Engineering rationale: a true rest-to-cruise transition requires near-
% saturated thrust + ~83 deg body rotation simultaneously, leaving zero
% margin for the cascade. After multiple iterations we adopt a fully
% rest-to-rest profile that is robustly trackable:
%
%   1. Vertical climb 0 -> 20 m            (rest-to-rest, pitch=90 deg)
%   2. Hover hold at 20 m                  (20 s)
%   3. Vertical climb 20 -> 1000 m         (rest-to-rest, pitch=90 deg)
%   4. Slow horizontal cruise NORTH        (rest-to-rest at 1000 m)
%
% The aircraft stays in tailsitter (hover) orientation throughout the climb
% phases. Horizontal motion at 1000 m is done as a slow rest-to-rest
% maneuver (NOT at V_LDmax — that would require sustained pitch transition
% which the cascade cannot perform robustly with current bandwidth).
%
% This profile sacrifices L/D-optimal cruise for robust trackability.
% Future work: implement L1 adaptive augmentation or NDI+INDI cascade
% redesign to enable true V_LDmax cruise.

if nargin < 2, opts = struct(); end
opts = set_default(opts, 'alt_low_m',        20);
opts = set_default(opts, 'alt_cruise_m',     1000);
opts = set_default(opts, 'cruise_distance_m', 1000);   % shorter, slow cruise
opts = set_default(opts, 'hover_duration_s',  20);
opts = set_default(opts, 'climb_safety',     2.0);    % 2x of physical minimum

g = 9.80665;
W = ac_cfg.mass * g;
rho = 1.225;

% --- Compute max L/D point (informational, NOT used as cruise target) ---
K = 1 / (pi * ac_cfg.wing.AR * ac_cfg.wing.e);
CL_LDmax = sqrt(ac_cfg.wing.CD0 / K);
V_LDmax = sqrt(2 * W / (rho * ac_cfg.wing.area * CL_LDmax));
alpha_LDmax = (CL_LDmax / ac_cfg.wing.CL_alpha) + ac_cfg.wing.alpha0;
LD_max = 0.5 / sqrt(ac_cfg.wing.CD0 * K);

info.V_LDmax = V_LDmax;
info.alpha_LDmax_deg = rad2deg(alpha_LDmax);
info.LD_max = LD_max;
info.CL_LDmax = CL_LDmax;

% --- Time allocation: use 30% of vertical envelope for high margin ---
%   a_envelope = N * T_max / m
%   a_vert_max = 0.3 * (a_envelope - g)   [conservative]
a_envelope = (ac_cfg.n_rotors * ac_cfg.rotor.thrust_max) / ac_cfg.mass;
a_vert_max = max(0.3 * (a_envelope - g), 0.2);   % minimum 0.2 m/s^2
peak_acc_coef = 15;

% Phase 1: vertical climb 0 -> -20 m
T_climb1 = max(8.0, opts.climb_safety * sqrt(peak_acc_coef * opts.alt_low_m / a_vert_max));

% Phase 2: hover hold at -20 m
T_hover = opts.hover_duration_s;

% Phase 3: vertical climb 20 -> 1000 m (rest-to-rest, slow)
d_alt = opts.alt_cruise_m - opts.alt_low_m;
T_climb2 = max(45.0, opts.climb_safety * sqrt(peak_acc_coef * d_alt / a_vert_max));

% Phase 4: rest-to-rest horizontal motion at 1000 m altitude
% Use horizontal envelope: a_horiz_max = 0.3 * sqrt(a_env^2 - g^2)
a_horiz_max = max(0.3 * sqrt(max(a_envelope^2 - g^2, 1.0)), 0.5);
T_cruise = max(30.0, opts.climb_safety * sqrt(peak_acc_coef * opts.cruise_distance_m / a_horiz_max));

% --- Build phases (rest-to-rest only, no rest-to-cruise) ---
p_origin   = [0; 0; 0];
p_alt_low  = [0; 0; -opts.alt_low_m];
p_alt_high = [0; 0; -opts.alt_cruise_m];
p_north    = [opts.cruise_distance_m; 0; -opts.alt_cruise_m];

ph1 = MissionTrajectory.make_rest_to_rest(p_origin,   p_alt_low,  T_climb1);
ph2 = MissionTrajectory.make_hover(p_alt_low, T_hover);
ph3 = MissionTrajectory.make_rest_to_rest(p_alt_low,  p_alt_high, T_climb2);
ph4 = MissionTrajectory.make_rest_to_rest(p_alt_high, p_north,    T_cruise);

traj = MissionTrajectory({ph1, ph2, ph3, ph4});

info.T_climb1  = T_climb1;
info.T_hover   = T_hover;
info.T_climb2  = T_climb2;
info.T_cruise  = T_cruise;
info.T_total   = traj.total_time();
info.distance_north_total = opts.cruise_distance_m;

% --- Print mission summary ---
fprintf('\n========================================================\n');
fprintf('   Conservative Realistic Mission Profile\n');
fprintf('========================================================\n');
fprintf('   Design point reference (informational, not used as target):\n');
fprintf('     V_LDmax     = %.2f m/s (%.0f km/h)\n', V_LDmax, V_LDmax*3.6);
fprintf('     alpha_LDmax = %.2f deg\n', rad2deg(alpha_LDmax));
fprintf('     L/D_max     = %.2f\n', LD_max);
fprintf('\n');
fprintf('   Trajectory (rest-to-rest only, robust tracking):\n');
fprintf('     a_vert_max  = %.2f m/s^2 (30%% of envelope)\n', a_vert_max);
fprintf('     a_horiz_max = %.2f m/s^2 (30%% of envelope)\n', a_horiz_max);
fprintf('\n');
fprintf('   Phase durations:\n');
fprintf('     1. Vertical climb 0->%.0fm    : %.1f s\n', opts.alt_low_m, T_climb1);
fprintf('     2. Hover at %.0f m              : %.1f s\n', opts.alt_low_m, T_hover);
fprintf('     3. Vertical climb %.0f->%.0fm  : %.1f s\n', opts.alt_low_m, opts.alt_cruise_m, T_climb2);
fprintf('     4. Slow cruise %.0fm north      : %.1f s\n', opts.cruise_distance_m, T_cruise);
fprintf('     TOTAL                          : %.1f s (%.2f min)\n', ...
    info.T_total, info.T_total/60);
fprintf('========================================================\n\n');

end


function s = set_default(s, field, default_val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = default_val;
    end
end
