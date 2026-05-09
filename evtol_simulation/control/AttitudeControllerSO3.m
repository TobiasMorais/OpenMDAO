classdef AttitudeControllerSO3 < handle
% ATTITUDECONTROLLERSO3  Geometric attitude tracking on SO(3).
%
% Reference: Lee, Leok, McClamroch (2010) "Geometric tracking control of a
% quadrotor UAV on SE(3)". Extended with integral term on so(3) for
% disturbance rejection (Goodarzi 2014).
%
% Error signals (Lie-algebra projections):
%   e_R  = 0.5 * vee( Rd^T R - R^T Rd )      attitude error (in body frame)
%   e_W  = omega - R^T Rd Wd                 angular-rate error
%   ei  += dt * e_R                           integral on so(3)
%
% Control torque (body frame):
%   M = -kR * e_R - kW * e_W - kI * ei
%       + omega x J omega
%       - J ( hat(omega) R^T Rd Wd  -  R^T Rd Wd_dot )
% (the feed-forward terms make the closed loop track Rd(t) globally.)
%
% Singularity-free, no Euler conversion, no quaternion unwinding.

    properties
        cfg
        ei = zeros(3,1);
        prev_t = NaN;
    end

    methods
        function obj = AttitudeControllerSO3(cfg)
            obj.cfg = cfg;
        end

        function reset(obj)
            obj.ei(:) = 0;
            obj.prev_t = NaN;
        end

        function M_cmd = compute(obj, R, omega, Rd, Wd, Wd_dot, J, t)
            % R: current body->world rotation
            % omega: body angular rates
            % Rd, Wd, Wd_dot: desired attitude, rate, rate_dot
            % J: inertia tensor
            e_R = so3_utils('errMat', R, Rd);

            Wd_b = R' * Rd * Wd;     % desired angular rate expressed in body frame
            e_W  = omega - Wd_b;

            % Integral term with anti-windup
            if isnan(obj.prev_t)
                dt = 0;
            else
                dt = max(0, t - obj.prev_t);
            end
            obj.prev_t = t;
            obj.ei = obj.ei + dt * e_R;
            obj.ei = max(-obj.cfg.so3.I_max, min(obj.cfg.so3.I_max, obj.ei));

            % Feed-forward terms
            term_ff = - J * ( so3_utils('hat', omega) * R' * Rd * Wd ...
                              - R' * Rd * Wd_dot );
            M_cmd =  - obj.cfg.so3.kR  * e_R ...
                     - obj.cfg.so3.kOm * e_W ...
                     - obj.cfg.so3.kI  * obj.ei ...
                     + cross(omega, J*omega) ...
                     + term_ff;
        end
    end
end
