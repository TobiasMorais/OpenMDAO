classdef PositionControllerNMPC < handle
% POSITIONCONTROLLERNMPC  Nonlinear MPC for outer position/velocity loop.
%
% Decision variable: sequence of body-frame specific-force commands a_cmd[k]
% over horizon N. Predictive model is point-mass under gravity:
%
%   p[k+1] = p[k] + dt * v[k]
%   v[k+1] = v[k] + dt * (R(q_des[k]) * a_cmd[k] + g_NED)
%
% where the desired thrust vector points from current state toward reference
% trajectory. Actually we work directly with NED-frame specific-force vector
% f_cmd, decoupling outer optimization from attitude:
%
%   v[k+1] = v[k] + dt * (f_cmd[k] + g_NED)
%
% Cost:
%   J = sum_k ( (p-p_ref)' Qp (p-p_ref) + (v-v_ref)' Qv (v-v_ref) + f' R f )
%       + e_N' Qf e_N
%
% Constraints:
%   |f_cmd[k]| <= a_max
%   tilt of f_cmd[k] from -Z_NED <= tilt_max  (only enforced once descent is needed)
%
% Outputs to inner loop:
%   F_cmd_NED (first step f_cmd[0]) -> converted to (T_total, qd) by mapper.
% This separation is the clean cascade GNC pattern.

    properties
        cfg
        last_solution = [];
        m_total
        % Disturbance estimator: integral of position error gives steady-state
        % offset compensation for unmodeled forces (aero drag, slipstream, wind).
        d_hat = zeros(3,1);     % NED-frame disturbance specific-force estimate
        d_max = 8.0;             % anti-windup limit [m/s^2]
    end

    methods
        function obj = PositionControllerNMPC(cfg, mass)
            obj.cfg = cfg;
            obj.m_total = mass;
        end

        function reset(obj)
            obj.last_solution = [];
            obj.d_hat = zeros(3,1);
        end

        function f_cmd = compute(obj, p, v, p_ref_traj, v_ref_traj)
            % p_ref_traj, v_ref_traj: 3 x (N+1) reference horizon
            N  = obj.cfg.nmpc.N;
            dt = obj.cfg.nmpc.dt;

            % --- Disturbance update (integrating action on position error) ---
            % d_hat += k_d * (p_ref - p) * dt_outer
            % This compensates steady-state offsets from aero/slipstream/wind.
            err_p = p_ref_traj(:,1) - p;
            k_d = 0.4;                    % integration gain
            obj.d_hat = obj.d_hat + k_d * err_p * (1/obj.cfg.f_outer);
            obj.d_hat = max(-obj.d_max, min(obj.d_max, obj.d_hat));

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
            % Apply disturbance compensation outside the optimization horizon
            % (treated as a known additive offset in the actuator).
            f_cmd = u_opt(1:3) + obj.d_hat;

            % Saturate combined command at a_max
            mag = norm(f_cmd);
            if mag > a_max
                f_cmd = f_cmd * (a_max / mag);
            end
        end

        function J = cost_and_dynamics(obj, u, p, v, pr, vr, N, dt)
            % Single-shooting cost evaluation
            J = 0;
            g_NED = [0;0;9.80665];
            for k = 1:N
                fk = u(3*(k-1)+1 : 3*k);
                v_next = v + dt * (fk + g_NED);
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
