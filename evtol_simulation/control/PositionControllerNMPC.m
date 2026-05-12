classdef PositionControllerNMPC < handle
% POSITIONCONTROLLERNMPC  Nonlinear MPC for outer position/velocity loop with
%                         offset-free disturbance observer.
%
% Decision variable: sequence of NED-frame specific-force commands f_cmd[k]
% over horizon N. Predictive model (offset-free):
%
%   p[k+1] = p[k] + dt * v[k]
%   v[k+1] = v[k] + dt * (f_cmd[k] + d_hat + g_NED)
%
% where d_hat is an integral disturbance estimate that captures unmodeled
% forces (aerodynamic drag, lift, slipstream download, wind, etc.). This is
% the standard "offset-free MPC" technique (Pannocchia & Rawlings 2003,
% Maeder & Morari 2010) — without it the controller commands trim assuming
% f = -g_NED but reality has f_real = -g_NED + d_real, producing a steady
% drift that the cascade saturates trying to compensate.
%
% Disturbance observer (closed-loop, anti-windup):
%   v_pred[k+1] = v[k] + dt * (f_cmd_applied[k] + d_hat[k] + g_NED)
%   error_v = v_actual[k+1] - v_pred[k+1]              (one-step prediction err)
%   d_hat[k+1] = d_hat[k] + k_obs * (error_v / dt)     (slow integration)
%   d_hat saturated to +/- d_max
%   Conditional integration: only update d_hat when actuator NOT saturated
%   in same direction (anti-windup).
%
% Cost: same as before, J = pos_err + vel_err + control effort.

    properties
        cfg
        last_solution = [];
        m_total
        % --- Disturbance observer state ---
        d_hat = zeros(3,1);     % NED-frame disturbance specific-force estimate
        d_max = 5.0;             % saturation [m/s^2]
        prev_v = [];             % velocity at previous outer-loop call
        prev_f_cmd = [];         % f_cmd commanded at previous call
        prev_t = NaN;            % time of previous call
        k_obs = 0.15;            % observer gain (slow integration)
    end

    methods
        function obj = PositionControllerNMPC(cfg, mass)
            obj.cfg = cfg;
            obj.m_total = mass;
        end

        function reset(obj)
            obj.last_solution = [];
            obj.d_hat = zeros(3,1);
            obj.prev_v = [];
            obj.prev_f_cmd = [];
            obj.prev_t = NaN;
        end

        function f_cmd = compute(obj, p, v, p_ref_traj, v_ref_traj)
            % p_ref_traj, v_ref_traj: 3 x (N+1) reference horizon
            N  = obj.cfg.nmpc.N;
            dt = obj.cfg.nmpc.dt;
            dt_outer = 1 / obj.cfg.f_outer;

            % --- Disturbance observer update ---
            % Compare actual velocity with one-step prediction from prev call
            if ~isempty(obj.prev_v) && ~isempty(obj.prev_f_cmd)
                v_pred = obj.prev_v + dt_outer * (obj.prev_f_cmd + obj.d_hat + [0;0;9.80665]);
                err_v = v - v_pred;
                % Integral update with conditional anti-windup:
                % don't increase |d_hat| if it's already saturated in same direction
                update = obj.k_obs * err_v;
                for ax = 1:3
                    if abs(obj.d_hat(ax)) >= obj.d_max
                        if sign(update(ax)) == sign(obj.d_hat(ax))
                            update(ax) = 0;   % saturated, don't push further
                        end
                    end
                end
                obj.d_hat = obj.d_hat + update;
                obj.d_hat = max(-obj.d_max, min(obj.d_max, obj.d_hat));
            end

            % Warm start
            if isempty(obj.last_solution)
                u0 = repmat([0;0;-9.80665], N, 1);    % hover hold
            else
                u0 = [obj.last_solution(4:end); obj.last_solution(end-2:end)];
            end

            % Bounds: |f| <= a_max along each axis is conservative box
            a_max = obj.cfg.nmpc.acc_max;
            lb = -a_max * ones(3*N, 1);
            ub =  a_max * ones(3*N, 1);

            cost_fn = @(u) obj.cost_and_dynamics(u, p, v, p_ref_traj, v_ref_traj, N, dt);

            try
                u_opt = fmincon(cost_fn, u0, [],[],[],[], lb, ub, ...
                                @(u) obj.tilt_constraint(u, N), obj.cfg.nmpc.opts);
            catch ME
                warning('NMPC fmincon failed (%s), using fallback PID', ME.message);
                u_opt = obj.fallback_pid(p, v, p_ref_traj, v_ref_traj, N);
            end

            obj.last_solution = u_opt;
            % f_cmd is the actuator command. d_hat is in the predictor (not added).
            f_cmd = u_opt(1:3);

            % Saturate at a_max (safety)
            mag = norm(f_cmd);
            if mag > a_max
                f_cmd = f_cmd * (a_max / mag);
            end

            % Save state for next call's observer update
            obj.prev_v = v;
            obj.prev_f_cmd = f_cmd;
        end

        function J = cost_and_dynamics(obj, u, p, v, pr, vr, N, dt)
            % Single-shooting cost evaluation with offset-free predictor.
            % d_hat is included so the planner accounts for unmodeled forces.
            J = 0;
            g_NED = [0;0;9.80665];
            for k = 1:N
                fk = u(3*(k-1)+1 : 3*k);
                v_next = v + dt * (fk + obj.d_hat + g_NED);
                p_next = p + dt * v;

                ep = p_next - pr(:, k+1);
                ev = v_next - vr(:, k+1);

                J = J + ep' * obj.cfg.nmpc.Q_pos * ep ...
                      + ev' * obj.cfg.nmpc.Q_vel * ev ...
                      + fk' * obj.cfg.nmpc.R_acc * fk;

                p = p_next; v = v_next;
            end
            % Terminal
            ef = [p - pr(:, N+1); v - vr(:, N+1)];
            J  = J + ef' * obj.cfg.nmpc.Qf * ef;
        end

        function [c, ceq] = tilt_constraint(obj, u, N)
            % Tilt of specific-force from local-vertical <= tilt_max
            tilt_max = obj.cfg.nmpc.tilt_max;
            c = zeros(N,1);
            for k = 1:N
                fk = u(3*(k-1)+1 : 3*k);
                fz = -fk(3);   % Z down: thrust pulls up => fz > 0 to oppose gravity
                fxy = norm(fk(1:2));
                c(k) = fxy - tan(tilt_max) * max(fz, 1.0);
            end
            ceq = [];
        end

        function u = fallback_pid(obj, p, v, pr, vr, N)
            % Simple PD if NMPC fails
            ep = pr(:, 2) - p;
            ev = vr(:, 2) - v;
            f0 = 1.5 * ep + 0.8 * ev - [0;0;9.80665];
            u = repmat(f0, N, 1);
        end
    end
end
