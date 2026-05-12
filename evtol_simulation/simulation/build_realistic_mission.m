function [traj, info] = build_realistic_mission(ac_cfg, opts)
% BUILD_REALISTIC_MISSION  Mission using DifferentialFlatness (proven in hover test).
%
% Engineering insight from many iterations: hover_only test PASSES
% with DifferentialFlatness, but ANY mission with MissionTrajectory
% diverges, even simple rest-to-rest climbs. To isolate the source of
% divergence, this version uses DifferentialFlatness with explicit
% waypoints — the same trajectory engine that powers the passing
% hover test.
%
% Mission profile (user request, conservatively timed):
%   1. Takeoff from -1 m to -20 m (climb 19 m)
%   2. Hover at -20 m for 20 s (same waypoint, long time_alloc)
%   3. Climb to -1000 m (slow)
%   4. Slow horizontal motion to 1000 m north at 1000 m altitude
%
% Initial state: aircraft starts at -1 m altitude (init_mode='hover_at_altitude'
% modified). The 1 m gap above ground avoids the motor-spin-up transient
% creating below-ground transient states that may confuse the cascade.

if nargin < 2, opts = struct(); end
opts = set_default(opts, 'alt_init_m',       1);
opts = set_default(opts, 'alt_low_m',        20);
opts = set_default(opts, 'alt_cruise_m',     1000);
opts = set_default(opts, 'cruise_distance_m', 1000);
opts = set_default(opts, 'hover_duration_s',  20);

g = 9.80665;
W = ac_cfg.mass * g;
rho = 1.225;

% --- Compute max L/D point (informational) ---
K = 1 / (pi * ac_cfg.wing.AR * ac_cfg.wing.e);
CL_LDmax = sqrt(ac_cfg.wing.CD0 / K);
V_LDmax = sqrt(2 * W / (rho * ac_cfg.wing.area * CL_LDmax));
alpha_LDmax = (CL_LDmax / ac_cfg.wing.CL_alpha) + ac_cfg.wing.alpha0;
LD_max = 0.5 / sqrt(ac_cfg.wing.CD0 * K);

info.V_LDmax = V_LDmax;
info.alpha_LDmax_deg = rad2deg(alpha_LDmax);
info.LD_max = LD_max;
info.CL_LDmax = CL_LDmax;

% --- Time allocation: 30% of envelope, 2x safety factor ---
a_envelope = (ac_cfg.n_rotors * ac_cfg.rotor.thrust_max) / ac_cfg.mass;
a_vert_max = max(0.3 * (a_envelope - g), 0.2);
a_horiz_max = max(0.3 * sqrt(max(a_envelope^2 - g^2, 1.0)), 0.5);
peak_acc_coef = 15;

T_climb1 = max(8.0, 2 * sqrt(peak_acc_coef * (opts.alt_low_m - opts.alt_init_m) / a_vert_max));
T_hover  = opts.hover_duration_s;
d_alt = opts.alt_cruise_m - opts.alt_low_m;
T_climb2 = max(45.0, 2 * sqrt(peak_acc_coef * d_alt / a_vert_max));
T_cruise = max(30.0, 2 * sqrt(peak_acc_coef * opts.cruise_distance_m / a_horiz_max));

% --- Waypoints (4 segments via 5 waypoints) ---
% NED: x=North, y=East, z=Down. psi=heading.
waypoints = [
    0, 0, -opts.alt_init_m,   0;    % start above ground
    0, 0, -opts.alt_low_m,    0;    % after climb to 20 m
    0, 0, -opts.alt_low_m,    0;    % after hover (same pt, long time)
    0, 0, -opts.alt_cruise_m, 0;    % after climb to 1000 m
    opts.cruise_distance_m, 0, -opts.alt_cruise_m, 0   % after slow cruise N
];
time_alloc = [T_climb1; T_hover; T_climb2; T_cruise];

% Use DifferentialFlatness (proven engine — same as hover_only test)
traj = DifferentialFlatness(waypoints, time_alloc);

info.T_climb1 = T_climb1;
info.T_hover  = T_hover;
info.T_climb2 = T_climb2;
info.T_cruise = T_cruise;
info.T_total  = traj.total_time();
info.distance_north_total = opts.cruise_distance_m;
info.init_mode = 'hover_at_altitude_init';
info.alt_init_m = opts.alt_init_m;

% --- Print summary ---
fprintf('\n========================================================\n');
fprintf('   Realistic Mission (DifferentialFlatness engine)\n');
fprintf('========================================================\n');
fprintf('   Design point reference:\n');
fprintf('     V_LDmax = %.1f m/s | alpha_LDmax = %.1f deg | L/D = %.1f\n', ...
    V_LDmax, rad2deg(alpha_LDmax), LD_max);
fprintf('\n');
fprintf('   Conservative trajectory (30%% envelope, 2x safety, all rest-to-rest):\n');
fprintf('     a_vert_max = %.2f m/s^2  |  a_horiz_max = %.2f m/s^2\n', ...
    a_vert_max, a_horiz_max);
fprintf('\n');
fprintf('   Initial state: aircraft at -%.0f m altitude (clear of ground)\n', opts.alt_init_m);
fprintf('   Phase durations:\n');
fprintf('     1. Climb %d->%dm     : %.1f s\n', opts.alt_init_m, opts.alt_low_m, T_climb1);
fprintf('     2. Hover at %d m       : %.1f s\n', opts.alt_low_m, T_hover);
fprintf('     3. Climb %d->%dm    : %.1f s\n', opts.alt_low_m, opts.alt_cruise_m, T_climb2);
fprintf('     4. Cruise %dm N      : %.1f s\n', opts.cruise_distance_m, T_cruise);
fprintf('     TOTAL                  : %.1f s (%.2f min)\n', info.T_total, info.T_total/60);
fprintf('========================================================\n\n');

end


function s = set_default(s, field, default_val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = default_val;
    end
end
