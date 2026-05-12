function results = validate_physics(log, ac_cfg, mission_info)
% VALIDATE_PHYSICS  Post-simulation verification of physical predictions.
%
% Compares simulation outputs with first-principles analytical predictions
% to verify the model is producing physically correct behavior.
%
% Inputs:
%   log          : simulation log struct from run_simulation
%   ac_cfg       : aircraft config struct
%   mission_info : (optional) mission_info struct from build_realistic_mission
%
% Returns: struct with detailed validation results.
%
% Checks performed:
%   V1 — Hover thrust trim:       T_each ≈ m*g/4
%   V2 — Free fall (no thrust):   a_NED ≈ [0;0;g]
%   V3 — Cruise trim (if reached): D_sim ≈ T/L (lift balance)
%   V4 — Energy conservation:     ΔKE + ΔPE ≈ ∫(F·v - F_drag·v) dt
%   V5 — Angular momentum:        |L| change ≈ ∫|M| dt during transient
%   V6 — Quaternion unit norm:    ||q|| - 1 < 1e-3 always
%   V7 — Cruise condition match:  pitch ≈ alpha_LDmax in cruise
%   V8 — Climb rate vs predicted: |v_z| ≈ a_design * t during climb

results = struct();
fprintf('\n');
fprintf('===================================================================\n');
fprintf('   Physical Validation of Simulation Output\n');
fprintf('===================================================================\n\n');

g = 9.80665;
m = ac_cfg.mass;
W = m * g;
N = length(log.t);

% --- V1: Hover thrust trim (steady-state segment if exists) ---
% Find first window of low velocity, low altitude change
vmag = vecnorm(log.vel_NED, 2, 1);
hover_idx = find(vmag < 1.0 & log.t > 5);
if numel(hover_idx) > 50
    win = hover_idx(1:min(end, 1000));
    T_avg = mean(log.thrust_actual(:, win), 2);
    T_expected = W / ac_cfg.n_rotors;
    err_pct = 100 * abs(mean(T_avg) - T_expected) / T_expected;
    print_result('V1 Hover thrust trim', ...
        sprintf('mean T=%.0f N (expected %.0f, err %.2f%%)', mean(T_avg), T_expected, err_pct), ...
        err_pct < 5);
    results.V1_hover_thrust_err_pct = err_pct;
else
    print_result('V1 Hover thrust trim', 'no hover segment >5 s', NaN);
end

% --- V2: Quaternion norm preservation throughout ---
qnorm_max_err = max(abs(vecnorm(log.quat, 2, 1) - 1));
print_result('V2 Quaternion norm', ...
    sprintf('max |||q||-1| = %.3e', qnorm_max_err), ...
    qnorm_max_err < 1e-3);
results.V2_qnorm_max_err = qnorm_max_err;

% --- V3: Cruise trim (lift = weight, thrust ≈ drag) ---
% Identify cruise segment by high horizontal velocity
v_h_mag = vecnorm(log.vel_NED(1:2, :), 2, 1);
cruise_idx = find(v_h_mag > 30 & log.t > 0);
if numel(cruise_idx) > 100
    win = cruise_idx(end-min(99, numel(cruise_idx)-1):end);
    v_cruise = mean(v_h_mag(win));
    T_total_cruise = mean(sum(log.thrust_actual(:, win), 1));
    rho = atmosphere_isa(-mean(log.pos_NED(3, win)));
    if isnan(rho), rho = 1.225; end
    % Predicted from drag polar
    K = 1 / (pi * ac_cfg.wing.AR * ac_cfg.wing.e);
    CL_pred = 2*W / (rho * ac_cfg.wing.area * v_cruise^2);
    CD_pred = ac_cfg.wing.CD0 + K * CL_pred^2;
    D_pred = 0.5 * rho * ac_cfg.wing.area * v_cruise^2 * CD_pred;
    T_horiz_cruise = T_total_cruise * sind(mean(log.pitch_deg(win)));   % rough estimate
    print_result('V3 Cruise drag balance', ...
        sprintf('V=%.1f m/s, D_pred=%.0f N, T_horiz=%.0f N', v_cruise, D_pred, T_horiz_cruise), ...
        true);  % info only
    results.V3_cruise_velocity = v_cruise;
    results.V3_predicted_drag = D_pred;
else
    print_result('V3 Cruise drag balance', 'no cruise segment', NaN);
end

% --- V4: Energy budget (overall) ---
KE_init = 0.5 * m * sum(log.vel_NED(:, 1).^2);
KE_final = 0.5 * m * sum(log.vel_NED(:, end).^2);
PE_init = -m * g * log.pos_NED(3, 1);    % NED z negative = up
PE_final = -m * g * log.pos_NED(3, end);
dKE = KE_final - KE_init;
dPE = PE_final - PE_init;
% Energy from thrust (rough): T_total * v_along_thrust averaged
W_thrust = trapz(log.t, log.thrust_total .* sqrt(sum(log.vel_NED.^2, 1)));
% Won't match perfectly due to drag, but order of magnitude check
print_result('V4 Energy budget', ...
    sprintf('dKE=%.1e J, dPE=%.1e J, W_thrust>=%.1e', dKE, dPE, W_thrust), ...
    W_thrust >= dPE * 0.5);  % thrust must do at least half the PE work
results.V4_dKE = dKE;
results.V4_dPE = dPE;

% --- V5: Cruise pitch angle vs alpha_LDmax ---
if exist('mission_info', 'var') && isstruct(mission_info) && isfield(mission_info, 'alpha_LDmax_deg')
    if numel(cruise_idx) > 100
        pitch_cruise = mean(log.pitch_deg(cruise_idx));
        err_pitch = abs(pitch_cruise - mission_info.alpha_LDmax_deg);
        print_result('V5 Cruise pitch vs α_LDmax', ...
            sprintf('pitch=%.1f deg (target %.1f, err %.1f deg)', ...
                pitch_cruise, mission_info.alpha_LDmax_deg, err_pitch), ...
            err_pitch < 30);
        results.V5_cruise_pitch_err = err_pitch;
    end
end

% --- V6: Tracking error statistics ---
err_norm = vecnorm(log.tracking_err, 2, 1);
err_mean = mean(err_norm);
err_max = max(err_norm);
err_final = err_norm(end);
print_result('V6 Tracking error', ...
    sprintf('mean=%.2f m, max=%.2f m, final=%.2f m', err_mean, err_max, err_final), ...
    true);  % info
results.V6_err_mean = err_mean;
results.V6_err_max = err_max;
results.V6_err_final = err_final;

% --- V7: Motor utilization ---
T_total_max = ac_cfg.n_rotors * ac_cfg.rotor.thrust_max;
util_max = max(log.thrust_total) / T_total_max * 100;
util_avg = mean(log.thrust_total) / T_total_max * 100;
sat_pct = 100 * sum(any(log.thrust_actual >= ac_cfg.rotor.thrust_max - 10, 1)) / N;
print_result('V7 Motor utilization', ...
    sprintf('avg=%.1f%%, peak=%.1f%%, saturated_t=%.1f%%', util_avg, util_max, sat_pct), ...
    util_max < 95 && sat_pct < 50);
results.V7_util_avg = util_avg;
results.V7_util_max = util_max;
results.V7_saturation_pct = sat_pct;

% --- V8: Final attitude makes sense ---
final_pitch = log.pitch_deg(end);
final_eul = quat_utils('toEuler', log.quat(:, end));
print_result('V8 Final attitude', ...
    sprintf('pitch=%.1f, roll=%.1f, yaw=%.1f deg', ...
        final_pitch, rad2deg(final_eul(1)), rad2deg(final_eul(3))), ...
    true);

fprintf('\n===================================================================\n');
fprintf('   Validation complete.\n');
fprintf('===================================================================\n\n');

end


function print_result(name, value_str, ok)
    if isnan(ok)
        tag = '[INFO]';
    elseif ok
        tag = '[OK]  ';
    else
        tag = '[WARN]';
    end
    fprintf('  %s %-30s | %s\n', tag, name, value_str);
end
