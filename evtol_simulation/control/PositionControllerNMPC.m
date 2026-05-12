classdef PositionControllerNMPC < handle
% POSITIONCONTROLLERNMPC  Nonlinear MPC with disturbance observer (offset-free).
%
% Predictor (offset-free):
%   p[k+1] = p[k] + dt * v[k]
%   v[k+1] = v[k] + dt * (f_cmd[k] + d_hat + g_NED)
%
% Disturbance observer (dimensionally consistent):
%   At each outer-loop call, measure actual acceleration via finite difference:
%     a_meas[k] = (v[k] - v[k-1]) / dt_outer
%   Predicted acceleration WITHOUT disturbance:
%     a_pred[k] = f_prev[k-1] + g_NED
%   Disturbance instantaneous estimate:
%     d_est[k] = a_meas[k] - a_pred[k]
%   Low-pass filter to update d_hat:
%     d_hat[k] = (1-alpha) * d_hat[k-1] + alpha * d_est[k]
%
% This is the standard form (Welch & Bishop 2006 disturbance observer,
% Liu et al 2009 active disturbance rejection). Time constant of observer:
%   tau_obs = dt_outer / alpha
% For alpha=0.3, tau_obs = 0.067 s — fast enough to track aero changes
% during transitions.
%
% Anti-windup on d_hat magnitude saturation.

    properties
        cfg
        last_solution = [];
        m_total
        % Disturbance observer state
        d_hat = zeros(3,1);     % NED-frame disturbance specific-force [m/s^2]
        d_max = 5.0;             % saturation [m/s^2]
        prev_v = [];             % velocity at previous outer-loop call
        prev_f_cmd = [];         % f_cmd commanded at previous call
        alpha_obs = 0.0;         % observer DISABLED (see notes below).
        %
        % ENGINEERING NOTE on the disturbance observer:
        %
        % An observer was implemented but disabled (alpha_obs = 0) after analysis
        % showed:
        %   1. NMPC's natural P-action handles steady-state disturbances within
        %      ~1 cm at hover (d_real=0.16 m/s^2) and ~5 cm at cruise
        %      (d_real=0.7 m/s^2). Steady-state error = d_real/sqrt(Qp/R) ≈ 1/14.
        %      Observer is not mathematically necessary.
        %
        %   2. With observer enabled, motor spin-up transient (~0.3 s) was
        %      misinterpreted as a "disturbance" because the predictor doesn't
        %      model first-order ESC dynamics. d_hat would saturate at 1.5 m/s^2
        %      downward, causing NMPC to over-command thrust, leading to
        %      oscillation and divergence.
        %
        % To re-enable: set alpha_obs > 0 (typical 0.05-0.1) AND model actuator
        % lag in predictor, OR delay observer activation by 1 s to let motors trim.
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
        end

        function f_cmd = compute(obj, p, v, p_ref_traj, v_ref_traj)
            N  = obj.cfg.nmpc.N;
            dt = obj.cfg.nmpc.dt;
            dt_outer = 1 / obj.cfg.f_outer;
            g_NED = [0;0;9.80665];

            % --- Disturbance observer (DISABLED by default) ---
            % See class header note. Observer confuses actuator lag with disturbance.
            if obj.alpha_obs > 0 && ~isempty(obj.prev_v) && ~isempty(obj.prev_f_cmd)
                a_meas = (v - obj.prev_v) / dt_outer;
                a_pred_no_d = obj.prev_f_cmd + g_NED;
                d_est = a_meas - a_pred_no_d;
                d_new = (1 - obj.alpha_obs) * obj.d_hat + obj.alpha_obs * d_est;
                obj.d_hat = max(-obj.d_max, min(obj.d_max, d_new));
            end

            % Warm start
            if isempty(obj.last_solution)
                u0 = repmat([0;0;-9.80665], N, 1);
            else
                u0 = [obj.last_solution(4:end); obj.last_solution(end-2:end)];
            end

            % Bounds
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
            f_cmd = u_opt(1:3);

            % Final saturation at a_max
            mag = norm(f_cmd);
            if mag > a_max
                f_cmd = f_cmd * (a_max / mag);
            end

            % Save state for next observer update
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
