classdef PositionControllerPDFF < handle
% POSITIONCONTROLLERPDFF  Classical PD controller with acceleration feedforward.
%
% This is the standard cascade outer loop for trajectory tracking of vectored-
% thrust UAVs (Mellinger & Kumar 2011, Bouabdallah 2007, Lee-Leok-McClamroch
% 2010). Used instead of NMPC because:
%   1. No optimization required — no fmincon convergence issues, no local minima
%   2. Robust to warm-start contamination (deterministic, stateless)
%   3. Handles moving references natively via the a_ref feedforward term
%   4. Well-understood frequency-domain properties and stability margins
%
% Control law (NED frame, specific force commanded):
%   f_cmd = -K_p (p - p_ref) - K_v (v - v_ref) + a_ref - g_NED
%
% With g_NED = [0;0;g], for hover at trim (p=p_ref, v=v_ref=0, a_ref=0):
%   f_cmd = -g_NED = [0;0;-g]   (specific thrust pointing UP in NED)
%
% Tuning: omega_n = 1 rad/s, zeta = 0.7. Bandwidth ratio with inner SO(3)
% (8 rad/s): 8:1 — solid cascade separation.
%
% Steady-state error from constant disturbance d:
%   err_ss = d / K_p = d / 1 = d  (in meters per m/s^2 of disturbance)
% For our aero drag d ~ 0.16 m/s^2, err_ss ~ 16 cm. Tighter K_p would reduce.

    properties
        K_p
        K_v
        a_max
        cfg
    end

    methods
        function obj = PositionControllerPDFF(cfg, mass)
            obj.cfg = cfg;
            % Tuning: critically-damped, omega_n = 1 rad/s
            % For stiffer tracking, increase. For less aggressive, decrease.
            omega_n = 1.0;
            zeta    = 0.85;
            obj.K_p = omega_n^2 * eye(3);                  % 1 * I
            obj.K_v = 2 * zeta * omega_n * eye(3);         % 1.7 * I
            obj.a_max = cfg.nmpc.acc_max;
        end

        function reset(obj)
            % Stateless controller
        end

        function f_cmd = compute(obj, p, v, p_ref_traj, v_ref_traj, a_ref)
            % p_ref_traj, v_ref_traj: 3 x (N+1) reference horizon
            % Only the FIRST column is used (current time reference).
            % a_ref: optional 3x1 acceleration feedforward (default 0).
            if nargin < 6 || isempty(a_ref)
                a_ref = zeros(3,1);
            end
            p_ref = p_ref_traj(:, 1);
            v_ref = v_ref_traj(:, 1);

            err_p = p - p_ref;
            err_v = v - v_ref;
            g_NED = [0; 0; 9.80665];

            f_cmd = -obj.K_p * err_p - obj.K_v * err_v + a_ref - g_NED;

            % Saturate at a_max (preserves direction)
            mag = norm(f_cmd);
            if mag > obj.a_max
                f_cmd = f_cmd * (obj.a_max / mag);
            end
        end
    end
end
