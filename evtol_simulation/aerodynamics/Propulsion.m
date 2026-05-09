classdef Propulsion < handle
% PROPULSION  Rotor model: BEMT-derived induced velocity + 1st-order ESC dynamics.
%
% --- Actuator dynamics (commanded vs actual rotor speed) ---
%   tau_m * \dot Omega_i + Omega_i = Omega_cmd_i
%
% --- Static thrust and torque (BEMT-consistent quadratic) ---
%   T_i = k_T * Omega_i^2
%   Q_i = k_Q * Omega_i^2     (sign-applied via spin direction for reaction torque)
%
% --- Induced velocity (Glauert/Froude momentum theory) ---
%   For axial inflow with free-stream component V_x along thrust axis:
%     v_i^2 + V_x * v_i - T / (2 * rho * A_disk) = 0
%   solution:
%     v_i = -V_x/2 + sqrt( (V_x/2)^2 + T/(2 rho A) )
%
% --- Slipstream velocity at downstream surface ---
%   V_slip = (k_s) * v_i      (k_s is wake development factor, ~1.5-2.0)
% Tube contraction sets the speed felt by surfaces immersed in the wake;
% added to local V_local in Aerodynamics.compute_forces_moments.
%
% --- Reaction torque on body ---
%   M_react = - sum_i sign(spin_i) * Q_i * \hat e_i
% (rotor accelerating CCW exerts CW torque on body, hence negative sign.)

    properties
        cfg
        rotor_axes          % 3 x N from Aircraft
        Omega_actual        % N x 1 current rotor speeds
    end

    methods
        function obj = Propulsion(cfg, rotor_axes)
            obj.cfg = cfg;
            obj.rotor_axes = rotor_axes;
            obj.Omega_actual = zeros(cfg.n_rotors, 1);
        end

        function step_actuator(obj, T_cmd, dt)
            % First-order ESC/motor dynamics on Omega.
            % T_cmd: 4x1 commanded thrust [N]
            T_cmd = max(0, min(obj.cfg.rotor.thrust_max, T_cmd));
            Omega_cmd = sqrt(T_cmd / obj.cfg.rotor.kT);
            tau = obj.cfg.rotor.tau_motor;
            obj.Omega_actual = obj.Omega_actual + (dt/tau) * (Omega_cmd - obj.Omega_actual);
            obj.Omega_actual = max(0, min(obj.cfg.rotor.omega_max, obj.Omega_actual));
        end

        function [F_prop_B, M_prop_B, v_i_per_rotor, Omega_signed] = compute(obj, V_B, rho)
            % Returns aggregate force/moment on body and per-rotor induced velocities.
            n = obj.cfg.n_rotors;
            F_prop_B = zeros(3,1);
            M_prop_B = zeros(3,1);
            v_i_per_rotor = zeros(n,1);
            Omega_signed  = zeros(n,1);

            for i = 1:n
                Om = obj.Omega_actual(i);
                T  = obj.cfg.rotor.kT * Om^2;
                Q  = obj.cfg.rotor.kQ * Om^2;

                axis_i = obj.rotor_axes(:, i);   % thrust direction in body
                r_i    = obj.cfg.rotor.position(i, :)';

                % Inflow component along thrust axis
                V_x = dot(V_B, axis_i);
                A   = obj.cfg.rotor.disk_area;

                if T > 1e-3
                    v_i = -V_x/2 + sqrt((V_x/2)^2 + T/(2*rho*A));
                else
                    v_i = 0;
                end
                v_i_per_rotor(i) = v_i;
                Omega_signed(i)  = obj.cfg.rotor.spin_dir(i) * Om;

                F_i = T * axis_i;
                M_thrust   = cross(r_i, F_i);
                M_reaction = -obj.cfg.rotor.spin_dir(i) * Q * axis_i;

                F_prop_B = F_prop_B + F_i;
                M_prop_B = M_prop_B + M_thrust + M_reaction;
            end
        end

        function V_slip_per_surf = slipstream_velocities(obj, v_i_per_rotor, V_B)
            % Distributes induced velocity into a 3 x N_surfaces matrix.
            % Surfaces order matches Aerodynamics.build_surface_strips().
            n_rotors = obj.cfg.n_rotors;
            n_surf   = size(obj.cfg.slipstream.bath_matrix, 2);
            V_slip_per_surf = zeros(3, n_surf);

            for s = 1:n_surf
                v_axial_total = 0;
                contributors  = 0;
                for i = 1:n_rotors
                    if obj.cfg.slipstream.bath_matrix(i, s) > 0
                        % contracted slipstream speed at downstream surface
                        v_axial_total = v_axial_total + obj.cfg.slipstream.k_s * v_i_per_rotor(i);
                        contributors  = contributors + 1;
                    end
                end
                if contributors > 0
                    v_avg = v_axial_total / contributors;
                    % Slipstream direction approximated as +X_B (rotor axes nominally +X_B)
                    % plus a fraction of free-stream so wing sees combined flow.
                    V_slip_per_surf(:, s) = [v_avg; 0; 0];
                end
            end
            % Note: This is a moderate-fidelity model. Swirl coupling and VLM mapping
            % can be plugged here by replacing the [v_avg;0;0] vector with a per-strip
            % swirl-corrected vector (see PDF section "Slipstream e Wing-Prop").
            %#ok<INUSD> V_B retained for future wake-skew extensions
        end

        function dT_max = thrust_rate_limit(obj)
            dT_max = obj.cfg.rotor.thrust_max / obj.cfg.rotor.tau_motor;
        end
    end
end
