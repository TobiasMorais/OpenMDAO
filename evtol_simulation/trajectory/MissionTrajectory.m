classdef MissionTrajectory < handle
% MISSIONTRAJECTORY  Composable multi-phase trajectory for realistic eVTOL missions.
%
% Supports four phase types:
%   1. 'rest_to_rest'  : 8th-order min-snap polynomial; v=a=j=0 at both ends
%   2. 'hover'         : constant position, zero velocity/accel
%   3. 'rest_to_cruise': 8th-order polynomial; v=0 at start, v=v_cruise at end
%   4. 'cruise'        : constant velocity, linear position
%   5. 'cruise_to_rest': 8th-order polynomial; v=v_cruise at start, v=0 at end
%
% Each phase carries its own (start_state, end_state, duration). Phase boundaries
% must satisfy continuity (next.start = prev.end), which the constructor verifies.
%
% Interface (drop-in replacement for DifferentialFlatness):
%   [pos, vel, acc, jerk, snap, psi_d, psi_dot] = traj.eval(t)
%   T = traj.total_time()

    properties
        phases       % cell array of phase structs
        T_total
        t_start_phase   % cumulative start time per phase
    end

    methods
        function obj = MissionTrajectory(phases)
            obj.phases = phases;
            % Compile coefficients for polynomial phases
            for k = 1:numel(phases)
                obj.phases{k} = obj.compile_phase(phases{k});
            end
            durations = cellfun(@(p) p.duration, obj.phases);
            obj.t_start_phase = [0; cumsum(durations(:))];
            obj.T_total = obj.t_start_phase(end);
            obj.verify_continuity();
        end

        function T = total_time(obj)
            T = obj.T_total;
        end

        function [pos, vel, acc, jerk, snap, psi_d, psi_dot] = eval(obj, t)
            t = max(0, min(t, obj.T_total));
            % Find phase
            k = find(t >= obj.t_start_phase(1:end-1) & t <= obj.t_start_phase(2:end), 1);
            if isempty(k), k = numel(obj.phases); end
            tau = t - obj.t_start_phase(k);

            ph = obj.phases{k};
            switch ph.type
                case 'hover'
                    pos = ph.position;
                    vel = zeros(3,1);
                    acc = zeros(3,1);
                    jerk = zeros(3,1);
                    snap = zeros(3,1);
                    psi_d = ph.psi;
                    psi_dot = 0;

                case 'cruise'
                    pos = ph.start_pos + ph.velocity * tau;
                    vel = ph.velocity;
                    acc = zeros(3,1);
                    jerk = zeros(3,1);
                    snap = zeros(3,1);
                    psi_d = ph.psi;
                    psi_dot = 0;

                case {'rest_to_rest', 'rest_to_cruise', 'cruise_to_rest'}
                    [pos, vel, acc, jerk, snap] = obj.eval_polynomial(ph, tau);
                    psi_d = ph.psi_start + (ph.psi_end - ph.psi_start) * (tau / ph.duration);
                    psi_dot = (ph.psi_end - ph.psi_start) / ph.duration;
            end
        end
    end

    % --- Static factory methods for common phases ---
    methods (Static)
        function ph = make_hover(position, duration, psi)
            if nargin < 3, psi = 0; end
            ph.type = 'hover';
            ph.position = position(:);
            ph.duration = duration;
            ph.psi = psi;
        end

        function ph = make_cruise(start_pos, velocity, duration, psi)
            if nargin < 4, psi = 0; end
            ph.type = 'cruise';
            ph.start_pos = start_pos(:);
            ph.velocity = velocity(:);
            ph.duration = duration;
            ph.end_pos = ph.start_pos + ph.velocity * duration;
            ph.psi = psi;
        end

        function ph = make_rest_to_rest(start_pos, end_pos, duration, psi_s, psi_e)
            if nargin < 4, psi_s = 0; end
            if nargin < 5, psi_e = psi_s; end
            ph.type = 'rest_to_rest';
            ph.start_pos = start_pos(:);
            ph.end_pos   = end_pos(:);
            ph.start_vel = zeros(3,1);
            ph.end_vel   = zeros(3,1);
            ph.start_acc = zeros(3,1);
            ph.end_acc   = zeros(3,1);
            ph.duration  = duration;
            ph.psi_start = psi_s;
            ph.psi_end   = psi_e;
        end

        function ph = make_rest_to_cruise(start_pos, end_pos, end_vel, duration, psi_s, psi_e)
            if nargin < 5, psi_s = 0; end
            if nargin < 6, psi_e = psi_s; end
            ph.type = 'rest_to_cruise';
            ph.start_pos = start_pos(:);
            ph.end_pos   = end_pos(:);
            ph.start_vel = zeros(3,1);
            ph.end_vel   = end_vel(:);
            ph.start_acc = zeros(3,1);
            ph.end_acc   = zeros(3,1);
            ph.duration  = duration;
            ph.psi_start = psi_s;
            ph.psi_end   = psi_e;
        end

        function ph = make_cruise_to_rest(start_pos, end_pos, start_vel, duration, psi_s, psi_e)
            if nargin < 5, psi_s = 0; end
            if nargin < 6, psi_e = psi_s; end
            ph.type = 'cruise_to_rest';
            ph.start_pos = start_pos(:);
            ph.end_pos   = end_pos(:);
            ph.start_vel = start_vel(:);
            ph.end_vel   = zeros(3,1);
            ph.start_acc = zeros(3,1);
            ph.end_acc   = zeros(3,1);
            ph.duration  = duration;
            ph.psi_start = psi_s;
            ph.psi_end   = psi_e;
        end
    end

    methods (Access = private)
        function ph = compile_phase(obj, ph)
            % Pre-compute polynomial coefficients for non-trivial phases.
            % 8 coefficients (8th-order polynomial) per dimension to satisfy
            % 8 boundary conditions: p, v, a, j at both ends.
            switch ph.type
                case {'rest_to_rest', 'rest_to_cruise', 'cruise_to_rest'}
                    T = ph.duration;
                    % Vandermonde-like system: rows = BCs at t=0 then t=T
                    A = [1 0 0 0 0 0 0 0;            % p(0)
                         0 1 0 0 0 0 0 0;            % p'(0)
                         0 0 2 0 0 0 0 0;            % p''(0)
                         0 0 0 6 0 0 0 0;            % p'''(0)
                         1 T T^2 T^3 T^4 T^5 T^6 T^7;
                         0 1 2*T 3*T^2 4*T^3 5*T^4 6*T^5 7*T^6;
                         0 0 2 6*T 12*T^2 20*T^3 30*T^4 42*T^5;
                         0 0 0 6 24*T 60*T^2 120*T^3 210*T^4];
                    coeffs = zeros(3, 8);
                    for d = 1:3
                        b = [ph.start_pos(d); ph.start_vel(d); ph.start_acc(d); 0;
                             ph.end_pos(d);   ph.end_vel(d);   ph.end_acc(d);   0];
                        coeffs(d, :) = (A \ b)';
                    end
                    ph.coeffs = coeffs;
            end
        end

        function [p, v, a, j, s] = eval_polynomial(obj, ph, tau)
            tau = max(0, min(tau, ph.duration));
            T = tau .^ (0:7)';
            Td = (1:7) .* tau .^ (0:6); Td = [0, Td]';
            Tdd = (1:6) .* (2:7) .* tau .^ (0:5); Tdd = [0, 0, Tdd]';
            Tddd = (1:5) .* (2:6) .* (3:7) .* tau .^ (0:4); Tddd = [0, 0, 0, Tddd]';
            Tdddd = (1:4) .* (2:5) .* (3:6) .* (4:7) .* tau .^ (0:3); Tdddd = [0, 0, 0, 0, Tdddd]';
            p = (ph.coeffs * T);
            v = (ph.coeffs * Td);
            a = (ph.coeffs * Tdd);
            j = (ph.coeffs * Tddd);
            s = (ph.coeffs * Tdddd);
        end

        function verify_continuity(obj)
            % Each phase's start must equal previous phase's end (position).
            for k = 2:numel(obj.phases)
                prev = obj.phases{k-1};
                curr = obj.phases{k};
                p_end_prev = obj.phase_end_position(prev);
                p_start_curr = obj.phase_start_position(curr);
                err = norm(p_end_prev - p_start_curr);
                if err > 1e-3
                    warning('MissionTrajectory: phase %d -> %d position discontinuity %.3f m', k-1, k, err);
                end
            end
        end

        function p = phase_end_position(obj, ph)
            switch ph.type
                case 'hover',  p = ph.position;
                case 'cruise', p = ph.end_pos;
                otherwise,     p = ph.end_pos;
            end
        end

        function p = phase_start_position(obj, ph)
            switch ph.type
                case 'hover',  p = ph.position;
                case 'cruise', p = ph.start_pos;
                otherwise,     p = ph.start_pos;
            end
        end
    end
end
