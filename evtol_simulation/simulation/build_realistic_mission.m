function [traj, info] = build_realistic_mission(ac_cfg, opts)
% BUILD_REALISTIC_MISSION  Realistic 5-phase tailsitter mission.
%
% Phases (NED frame, x=North, y=East, z=Down):
%
%   1. VERTICAL CLIMB to 20 m (rest-to-rest)
%      - Pure vertical from origin to [0;0;-20]
%      - Tailsitter pitch = +90 deg (body x = up) maintained
%      - Duration: bounded by 50% of vertical thrust envelope (a_vert <= ~1.1 m/s^2)
%
%   2. HOVER HOLD at 20 m for 20 s
%      - Position constant at [0;0;-20]
%      - Velocity = 0, acceleration = 0
%      - Tests stability of the cascade in steady-state
%
%   3. CLIMB-AND-TRANSITION to 1000 m altitude (rest-to-cruise)
%      - From [0;0;-20] @ rest to [d_climb;0;-1000] @ cruise velocity
%      - Aircraft transitions from hover (theta=90 deg) to cruise pitch
%      - Slow enough to respect physical pitch authority (~5 deg/s)
%      - Pre-position for cruise vector
%
%   4. STEADY CRUISE at 1000 m altitude, north direction
%      - Velocity = V_LDmax (max efficiency aerodynamic point)
%      - Pitch = alpha_LDmax (aircraft trims at max L/D AoA)
%      - 3 km north distance covered
%
%   5. (optional) CRUISE-TO-REST deceleration at 1000 m
%      - Bring vehicle back to hover for landing prep (out of scope here)
%
% Inputs:
%   ac_cfg : aircraft configuration struct
%   opts   : struct (optional) with fields:
%     .climb_rate_design  [m/s]   target vertical climb rate (default: 1.5)
%     .cruise_distance_m  [m]     horizontal cruise distance (default: 3000)
%     .hover_duration_s   [s]     hover hold duration (default: 20)
%     .alt_low_m          [m]     intermediate hover altitude (default: 20)
%     .alt_cruise_m       [m]     cruise altitude (default: 1000)
%
% Outputs:
%   traj : MissionTrajectory instance
%   info : struct with derived design values (V_LDmax, alpha_LDmax, etc.)

if nargin < 2, opts = struct(); end
opts = set_default(opts, 'climb_rate_design', 1.5);     % m/s vertical climb rate
opts = set_default(opts, 'cruise_distance_m', 3000);     % m
opts = set_default(opts, 'hover_duration_s', 20);        % s
opts = set_default(opts, 'alt_low_m', 20);               % m
opts = set_default(opts, 'alt_cruise_m', 1000);          % m

g = 9.80665;
W = ac_cfg.mass * g;
rho = 1.225;

% --- Compute max L/D cruise design point ---
K = 1 / (pi * ac_cfg.wing.AR * ac_cfg.wing.e);
CL_LDmax = sqrt(ac_cfg.wing.CD0 / K);
V_LDmax = sqrt(2 * W / (rho * ac_cfg.wing.area * CL_LDmax));
alpha_LDmax = (CL_LDmax / ac_cfg.wing.CL_alpha) + ac_cfg.wing.alpha0;
LD_max = 0.5 / sqrt(ac_cfg.wing.CD0 * K);

info.V_LDmax = V_LDmax;
info.alpha_LDmax_deg = rad2deg(alpha_LDmax);
info.LD_max = LD_max;
info.CL_LDmax = CL_LDmax;

% --- Phase 1: vertical climb 0 -> alt_low ---
% Conservative time using 50% of vertical envelope
T_per_max = 0.5 * ac_cfg.rotor.thrust_max;
a_envelope_full = (ac_cfg.n_rotors * ac_cfg.rotor.thrust_max) / ac_cfg.mass;
a_max_vert = max(0.5 * (a_envelope_full - g), 0.3);
peak_acc_coef = 15;
T_climb1 = max(8.0, sqrt(peak_acc_coef * opts.alt_low_m / a_max_vert));

ph1 = MissionTrajectory.make_rest_to_rest(...
    [0; 0; 0], [0; 0; -opts.alt_low_m], T_climb1);

% --- Phase 2: hover hold ---
ph2 = MissionTrajectory.make_hover([0; 0; -opts.alt_low_m], opts.hover_duration_s);

% --- Phase 3: climb-and-transition to cruise altitude ---
% End state: at cruise altitude, moving north at V_LDmax
% Position: must be far enough north that the transition is smooth
% We use a pre-positioning distance of (V_LDmax * T_transition / 2)
d_alt_to_cruise = opts.alt_cruise_m - opts.alt_low_m;
T_transition = max(45.0, sqrt(peak_acc_coef * d_alt_to_cruise / a_max_vert));
% While transitioning, accumulate horizontal distance:
% Approximate: average velocity = V_LDmax/2 over T_transition
d_north_during_transition = V_LDmax * T_transition / 2;

end_pos_3 = [d_north_during_transition; 0; -opts.alt_cruise_m];
end_vel_3 = [V_LDmax; 0; 0];
ph3 = MissionTrajectory.make_rest_to_cruise(...
    [0; 0; -opts.alt_low_m], end_pos_3, end_vel_3, T_transition);

% --- Phase 4: steady cruise at 1000 m ---
T_cruise = opts.cruise_distance_m / V_LDmax;
ph4 = MissionTrajectory.make_cruise(end_pos_3, end_vel_3, T_cruise);

% Build trajectory
traj = MissionTrajectory({ph1, ph2, ph3, ph4});

info.T_climb1     = T_climb1;
info.T_hover      = opts.hover_duration_s;
info.T_transition = T_transition;
info.T_cruise     = T_cruise;
info.T_total      = traj.total_time();
info.distance_north_total = end_pos_3(1) + opts.cruise_distance_m;

% --- Print mission summary ---
fprintf('\n========================================================\n');
fprintf('   Realistic Mission Profile\n');
fprintf('========================================================\n');
fprintf('   Aerodynamic design point (max L/D):\n');
fprintf('     V_LDmax     = %.2f m/s (%.0f km/h)\n', V_LDmax, V_LDmax*3.6);
fprintf('     alpha_LDmax = %.2f deg\n', rad2deg(alpha_LDmax));
fprintf('     CL_LDmax    = %.3f\n', CL_LDmax);
fprintf('     L/D_max     = %.2f\n', LD_max);
fprintf('   Phase durations:\n');
fprintf('     1. Vertical climb 0->%.0fm    : %.1f s\n', opts.alt_low_m, T_climb1);
fprintf('     2. Hover at %.0f m              : %.1f s\n', opts.alt_low_m, opts.hover_duration_s);
fprintf('     3. Climb+transition %d->%dm : %.1f s\n', ...
    opts.alt_low_m, opts.alt_cruise_m, T_transition);
fprintf('     4. Cruise %.0f km north         : %.1f s\n', opts.cruise_distance_m/1000, T_cruise);
fprintf('     TOTAL                         : %.1f s (%.2f min)\n', ...
    info.T_total, info.T_total/60);
fprintf('   Distance covered: %.0f m north, %.0f m altitude\n', ...
    info.distance_north_total, opts.alt_cruise_m);
fprintf('========================================================\n\n');

end


function s = set_default(s, field, default_val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = default_val;
    end
end
