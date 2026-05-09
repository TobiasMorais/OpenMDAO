classdef Aircraft < handle
% AIRCRAFT  Rigid-body 6-DOF tailsitter eVTOL with quaternion attitude.
%
% State vector x (13 components):
%   x(1:3)   = p_NED      [m]   position in inertial NED frame
%   x(4:6)   = v_NED      [m/s] velocity in NED
%   x(7:10)  = q          [-]   quaternion body->world (Hamilton, scalar first)
%   x(11:13) = omega_B    [rad/s] body angular rates [p; q; r]
%
% Newton-Euler equations of motion (body frame, generalized for J_xz coupling):
%
%   m * (\dot v^B + omega^B x v^B) = F_aero^B + F_prop^B + R_I^B * F_grav^I
%   J  *  \dot omega^B + omega^B x (J * omega^B) + M_gyro_rotor^B
%                                  = M_aero^B + M_prop^B
%
% Translation in NED frame is integrated as:
%   \dot p_NED = v_NED
%   \dot v_NED = (1/m) * [R(q) * (F_aero^B + F_prop^B)] + g_NED
%
% Quaternion kinematics (singularity-free):
%   \dot q = 0.5 * Omega(omega^B) * q
%
% Inertia matrix (PDF section "Modelagem de Massa..."):
%   J = [J_xx,    0,  -J_xz;
%          0,  J_yy,    0;
%       -J_xz,    0,  J_zz]
%
% Rotor gyroscopic torque (PDF section "Acoplamento Giroscópico"):
%   M_gyro = - omega_B x sum_i ( J_rotor_i * Omega_rotor_i * \hat e_i )
% where \hat e_i is unit thrust axis of rotor i in body frame (cant-rotated).

    properties
        cfg                 % aircraft config struct
        rotor_axes          % 3x4 matrix of unit thrust axes in body frame (cant-applied)
        gravity = 9.80665;
    end

    methods
        function obj = Aircraft(cfg)
            obj.cfg = cfg;
            obj.rotor_axes = obj.compute_rotor_axes();
        end

        function axes = compute_rotor_axes(obj)
            % Build unit thrust direction for each rotor by applying
            % yaw_cant about Z_B then pitch_cant about Y_B to nominal +X_B.
            n = obj.cfg.n_rotors;
            axes = zeros(3, n);
            for i = 1:n
                ay = obj.cfg.rotor.yaw_cant(i);
                ap = obj.cfg.rotor.pitch_cant(i);
                Rz = [cos(ay), -sin(ay), 0;
                      sin(ay),  cos(ay), 0;
                            0,        0, 1];
                Ry = [ cos(ap), 0, sin(ap);
                             0, 1,       0;
                      -sin(ap), 0, cos(ap)];
                axes(:, i) = Rz * Ry * [1; 0; 0];
            end
        end

        function xdot = dynamics(obj, x, F_aero_B, M_aero_B, F_prop_B, M_prop_B, M_gyro_B)
            % Integrated EOM. F/M arguments are pre-computed by Aerodynamics/Propulsion modules.
            v_NED = x(4:6);
            q     = quat_utils('norm', x(7:10));
            w_B   = x(11:13);

            R_BW = quat_utils('toR', q);          % body -> world
            R_WB = R_BW';                          % world -> body

            % Translation (NED frame integration is cleaner for navigation)
            F_total_B = F_aero_B + F_prop_B;
            a_NED = (R_BW * F_total_B) / obj.cfg.mass + [0; 0; obj.gravity];

            % Rotation (body frame)
            J = obj.cfg.J;
            M_total_B = M_aero_B + M_prop_B - M_gyro_B - cross(w_B, J * w_B);
            wdot_B    = obj.cfg.J_inv * M_total_B;

            % Quaternion kinematics
            qdot = quat_utils('kinematic', q, w_B);

            xdot = [v_NED; a_NED; qdot; wdot_B];
            % Note: R_WB unused here but available if a body-frame velocity formulation
            % is preferred (PDF uses both forms; we picked NED for navigation simplicity).
            %#ok<NASGU>
        end

        function M = compute_rotor_gyro(obj, omega_rotors_signed, w_B)
            % M_gyro = omega_B x H_rotor,    H_rotor = sum_i J_r * Omega_signed_i * \hat e_i
            % omega_rotors_signed: 4x1 with sign already applied (CCW=+, CW=-)
            J_r  = obj.cfg.rotor.J_rotor;
            H = zeros(3,1);
            for i = 1:obj.cfg.n_rotors
                H = H + J_r * omega_rotors_signed(i) * obj.rotor_axes(:, i);
            end
            M = cross(w_B, H);
        end

        function x0 = initial_state(obj, mode)
            % Build initial state vector for canonical scenarios.
            % mode: 'hover'   -> tailsitter pointing zenith (theta = +90 deg)
            %       'cruise'  -> level flight at cruise speed
            switch lower(mode)
                case 'hover'
                    p0 = [0; 0; 0];
                    v0 = [0; 0; 0];
                    q0 = quat_utils('fromEuler', 0, deg2rad(90), 0);
                    w0 = [0; 0; 0];
                case 'cruise'
                    p0 = [0; 0; -50];           % 50 m altitude (NED z down)
                    v0 = [50; 0; 0];            % 50 m/s North
                    q0 = quat_utils('fromEuler', 0, 0, 0);
                    w0 = [0; 0; 0];
                otherwise
                    error('Unknown mode: %s', mode);
            end
            x0 = [p0; v0; q0; w0];
        end
    end
end
