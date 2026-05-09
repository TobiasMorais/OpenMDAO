classdef AttitudeControllerINDI < handle
% ATTITUDECONTROLLERINDI  Incremental Nonlinear Dynamics Inversion (INDI).
%
% Reference: Smeur, Chu, de Croon (2016) "Adaptive Incremental Nonlinear Dynamics
% Inversion for Attitude Control of Micro Air Vehicles". Adapted for tailsitter.
%
% Core principle: instead of inverting the full nonlinear model, we use the
% latest filtered angular acceleration measurement and only invert the
% incremental relationship between control inputs and acceleration.
%
%   omega_dot = f(x) + G * u
%   omega_dot_filt   <-  low-pass filter of measured dot
%   u_filt           <-  low-pass filter of actuator state
%
%   Delta_u = G^+ * ( nu - omega_dot_filt )
%   u_cmd   = u_filt + Delta_u
%
% where nu is the desired angular acceleration produced by an outer
% attitude-tracking law (here: P on quaternion error, then P on rate error):
%
%   omega_des = K_att * vec(qe)
%   nu        = K_rate * (omega_des - omega) + omega_des_dot
%
% Robustness: G is a constant nominal effectiveness Jacobian (3 x m, with m
% = number of actuators); we use it for the incremental update only. As long as
% the measurement of omega_dot is accurate and the control loop is fast,
% modeling errors in f(x) cancel out.

    properties
        cfg
        G                   % effectiveness matrix d(omega_dot)/du, 3 x m
        omega_dot_filt = zeros(3,1);
        u_filt
        prev_omega = zeros(3,1);
        prev_t = NaN;
    end

    methods
        function obj = AttitudeControllerINDI(cfg, G_init, n_actuators)
            obj.cfg = cfg;
            obj.G = G_init;
            obj.u_filt = zeros(n_actuators, 1);
        end

        function reset(obj)
            obj.omega_dot_filt(:) = 0;
            obj.u_filt(:) = 0;
            obj.prev_omega(:) = 0;
            obj.prev_t = NaN;
        end

        function [Delta_u, nu] = compute(obj, q, omega, qd, omega_des_dot, u_meas, t)
            % q: current quaternion, qd: desired
            % omega: current body rates, omega_des_dot: feed-forward angular accel
            % u_meas: latest actuator state vector
            qe = quat_utils('errMul', q, qd);
            % Reduced attitude error (small-angle vector part)
            e_att  = qe(2:4);

            omega_des = - obj.cfg.indi.K_att * e_att;
            e_rate    = omega_des - omega;

            nu = obj.cfg.indi.K_rate * e_rate + omega_des_dot;

            % Update filtered omega_dot from finite difference + LP
            if ~isnan(obj.prev_t) && t > obj.prev_t
                dt = t - obj.prev_t;
                omega_dot_meas = (omega - obj.prev_omega) / dt;
                a = obj.cfg.indi.lp_omega * dt;
                a = a / (1 + a);
                obj.omega_dot_filt = (1-a)*obj.omega_dot_filt + a*omega_dot_meas;

                au = obj.cfg.indi.lp_u * dt; au = au/(1+au);
                obj.u_filt = (1-au)*obj.u_filt + au*u_meas;
            else
                obj.omega_dot_filt = zeros(3,1);
                obj.u_filt = u_meas;
            end
            obj.prev_omega = omega;
            obj.prev_t = t;

            % Incremental inversion via weighted pseudoinverse of G
            Gp = pinv(obj.G);
            Delta_u = Gp * (nu - obj.omega_dot_filt);
        end

        function update_effectiveness(obj, G_new)
            obj.G = G_new;
        end
    end
end
