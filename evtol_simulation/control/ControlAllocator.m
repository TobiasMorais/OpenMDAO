classdef ControlAllocator < handle
% CONTROLALLOCATOR  Maps virtual controls to actuator commands.
%
% Virtual controls (6-vector):
%   nu = [F_x_B; F_y_B; F_z_B; M_x_B; M_y_B; M_z_B]
%
% Actuators (7-vector, can extend if surfaces add):
%   u = [T_1; T_2; T_3; T_4; delta_e; delta_a; delta_r]
%
% Effectiveness matrix B (6x7) is built from:
%   - rotor thrust direction \hat e_i (with cant) and position r_i
%   - control surface eta * S_eff * arm contributions
%
% Force/moment per rotor i for unit thrust:
%     dF/dT_i = \hat e_i
%     dM/dT_i = r_i x \hat e_i
%
% Allocation via weighted pseudoinverse:
%     u = W^-1 B' (B W^-1 B')^-1 nu
% with box constraints projected via active-set.

    properties
        cfg
        ac_cfg
        rotor_axes
        B               % effectiveness matrix
    end

    methods
        function obj = ControlAllocator(ctrl_cfg, ac_cfg, rotor_axes)
            obj.cfg = ctrl_cfg;
            obj.ac_cfg = ac_cfg;
            obj.rotor_axes = rotor_axes;
            obj.B = obj.build_effectiveness();
        end

        function B = build_effectiveness(obj)
            % 6 x 7 effectiveness mapping
            B = zeros(6, 7);
            for i = 1:obj.ac_cfg.n_rotors
                e_i = obj.rotor_axes(:, i);
                r_i = obj.ac_cfg.rotor.position(i, :)';
                B(1:3, i) = e_i;                 % force per unit thrust
                B(4:6, i) = cross(r_i, e_i);     % moment per unit thrust
                % Reaction torque ignored at allocation level (small);
                % could be added: B(4:6, i) -= sign_i * (kQ/kT) * e_i  for full coupling.
            end

            % Surface effectiveness (linearized about cruise q_bar; allocator scales by speed).
            S = obj.ac_cfg.htail.area;
            arm = obj.ac_cfg.htail.arm;
            % Elevator: pitch moment (col 5)
            B(5, 5) = -obj.ac_cfg.htail.elev_eff * S * arm;
            % Aileron: roll moment (col 6)
            B(4, 6) = obj.ac_cfg.ail.eff * obj.ac_cfg.wing.area * obj.ac_cfg.wing.span;
            % Rudder: yaw moment (col 7)
            B(6, 7) = -obj.ac_cfg.vtail.rud_eff * obj.ac_cfg.vtail.area * obj.ac_cfg.vtail.arm;
        end

        function [u, info] = allocate(obj, nu, q_bar)
            % nu: 6x1 virtual command, q_bar: dynamic pressure (scales surfaces)
            B_eff = obj.B;
            % Scale aerodynamic columns by q_bar so surfaces produce real force/moment
            B_eff(:, 5:7) = B_eff(:, 5:7) * max(q_bar, 1.0);

            W = obj.cfg.alloc.W;
            % Weighted pseudoinverse: u = W^-1 B' (B W^-1 B')^-1 nu
            Wi = inv(W);
            BWB = B_eff * Wi * B_eff';
            BWB = BWB + 1e-3 * eye(size(BWB));   % regularize
            u_unc = Wi * B_eff' * (BWB \ nu);

            % Saturate
            u = u_unc;
            T_min = obj.cfg.alloc.T_min;
            T_max = obj.cfg.alloc.T_max;
            u(1:4) = max(T_min, min(T_max, u(1:4)));
            u(5)   = max(-obj.ac_cfg.surf.elev_max, min(obj.ac_cfg.surf.elev_max, u(5)));
            u(6)   = max(-obj.ac_cfg.ail.delta_max, min(obj.ac_cfg.ail.delta_max, u(6)));
            u(7)   = max(-obj.ac_cfg.surf.rud_max,  min(obj.ac_cfg.surf.rud_max,  u(7)));

            info.unconstrained = u_unc;
            info.residual = nu - B_eff * u;
            info.B = B_eff;
        end

        function G = attitude_effectiveness(obj)
            % 3 x 4: d(omega_dot)/dT_i for INDI (rotors only)
            % omega_dot = J^-1 (M_thrust) ; M_thrust(:,i) = r_i x e_i
            n = obj.ac_cfg.n_rotors;
            G = zeros(3, n);
            for i = 1:n
                M_per_T = cross(obj.ac_cfg.rotor.position(i,:)', obj.rotor_axes(:, i));
                G(:, i) = obj.ac_cfg.J_inv * M_per_T;
            end
        end
    end
end
