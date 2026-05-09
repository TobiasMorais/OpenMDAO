classdef DifferentialFlatness < handle
% DIFFERENTIALFLATNESS  Polynomial trajectory in flat outputs (x,y,z,psi).
%
% A vectored-thrust eVTOL is differentially flat in (sigma) = [x; y; z; psi].
% Once a smooth sigma(t) up to 4th derivative is given, the body attitude,
% thrust, and angular rates are recovered algebraically:
%
%   T (vec)  =  m * (sigma_ddot - g_NED)        (specific thrust in NED)
%   x_B      =  T / |T|
%   y_B      =  ... (heading constraint)
%   omega_B  =  recovered from third derivative (jerk)
%   omega_dot_B = recovered from fourth derivative (snap)
%
% Trajectory class: minimum-snap polynomials between waypoints
% (Mellinger & Kumar 2011). Each segment k between waypoints uses 8th-order
% polynomial in t to make all derivatives up to snap continuous.

    properties
        waypoints       % N x 4: [x y z psi] at each waypoint
        time_alloc      % (N-1) x 1 time per segment
        coeffs          % 4 x 8 x (N-1)  [x;y;z;psi] coefficients per segment
    end

    methods
        function obj = DifferentialFlatness(waypoints, time_alloc)
            obj.waypoints = waypoints;
            obj.time_alloc = time_alloc;
            obj.coeffs = obj.solve_minsnap();
        end

        function C = solve_minsnap(obj)
            % Solve QP for min-snap polynomials.
            % For brevity we use a closed-form 8th-order polynomial that
            % matches position, vel, accel, jerk = 0 at endpoints (rest-to-rest).
            N_seg = size(obj.waypoints, 1) - 1;
            C = zeros(4, 8, N_seg);

            for k = 1:N_seg
                T = obj.time_alloc(k);
                p0 = obj.waypoints(k,   :);
                p1 = obj.waypoints(k+1, :);
                for d = 1:4
                    % 8th-order polynomial: p(t) = sum_{j=0..7} c_j * t^j
                    % BCs: p(0)=p0, p'(0)=0, p''(0)=0, p'''(0)=0
                    %      p(T)=p1, p'(T)=0, p''(T)=0, p'''(T)=0
                    % Closed-form result yields a symmetric polynomial.
                    A = [1 0 0 0 0 0 0 0;
                         0 1 0 0 0 0 0 0;
                         0 0 2 0 0 0 0 0;
                         0 0 0 6 0 0 0 0;
                         1 T T^2 T^3 T^4 T^5 T^6 T^7;
                         0 1 2*T 3*T^2 4*T^3 5*T^4 6*T^5 7*T^6;
                         0 0 2 6*T 12*T^2 20*T^3 30*T^4 42*T^5;
                         0 0 0 6 24*T 60*T^2 120*T^3 210*T^4];
                    b = [p0(d); 0; 0; 0; p1(d); 0; 0; 0];
                    C(d, :, k) = (A \ b)';
                end
            end
        end

        function [pos, vel, acc, jerk, snap, psi_d, psi_dot] = eval(obj, t)
            % Evaluate trajectory at time t (0 <= t <= sum(time_alloc)).
            T_cumul = [0; cumsum(obj.time_alloc(:))];
            t = max(0, min(t, T_cumul(end)));
            k = find(t <= T_cumul, 1) - 1;
            if isempty(k) || k < 1, k = 1; end
            tau = t - T_cumul(k);

            powers   = tau .^ (0:7)';
            d1coef   = (0:7)';
            d2coef   = (0:7)' .* max(0, (-1:6)');
            d3coef   = (0:7)' .* max(0, (-1:6)') .* max(0, (-2:5)');
            d4coef   = (0:7)' .* max(0, (-1:6)') .* max(0, (-2:5)') .* max(0, (-3:4)');

            sigma   = zeros(4,1);
            sigmad  = zeros(4,1);
            sigmadd = zeros(4,1);
            sigmaddd  = zeros(4,1);
            sigmadddd = zeros(4,1);

            for d = 1:4
                c = squeeze(obj.coeffs(d, :, k))';
                sigma(d)   = c' * powers;
                sigmad(d)  = sum(d1coef .* c .* (tau.^max(0,(-1:6))'));
                sigmadd(d) = sum(d2coef .* c .* (tau.^max(0,(-2:5))'));
                sigmaddd(d) = sum(d3coef .* c .* (tau.^max(0,(-3:4))'));
                sigmadddd(d) = sum(d4coef .* c .* (tau.^max(0,(-4:3))'));
            end

            pos    = sigma(1:3);
            vel    = sigmad(1:3);
            acc    = sigmadd(1:3);
            jerk   = sigmaddd(1:3);
            snap   = sigmadddd(1:3);
            psi_d  = sigma(4);
            psi_dot = sigmad(4);
        end

        function [F_des_NED, qd, Wd] = recover_states(obj, t, mass, g)
            % Recover thrust vector + desired attitude from flat outputs.
            [~, ~, acc, jerk, ~, psi_d, ~] = obj.eval(t);
            F_des_NED = mass * (acc - [0;0;g]);   % g_NED has +Z down so subtract
            % qd via force_to_attitude logic (heading psi_d, thrust along F_des)
            % We replicate the mapping inline here for clarity.
            Fmag = norm(F_des_NED);
            if Fmag < 1e-3
                qd = quat_utils('fromEuler', 0, deg2rad(90), psi_d);
                Wd = zeros(3,1);
                return;
            end
            x_B = F_des_NED / Fmag;
            h = [cos(psi_d); sin(psi_d); 0];
            y_B = cross([0;0;-1], x_B);
            if norm(y_B) < 1e-3, y_B = cross(h, x_B); end
            y_B = y_B / norm(y_B);
            R = [x_B, y_B, cross(x_B, y_B)];
            qd = quat_utils('fromR', R);

            % Body angular velocity: derive from jerk projection (Mellinger 2011)
            % omega = R^T * (h_w x x_B) / |T|, where h_w = m * jerk
            h_w = mass * jerk;
            num = cross(h_w, x_B);
            Wd_world = num / max(Fmag, 1e-3);
            Wd = R' * Wd_world;
        end

        function T_total = total_time(obj)
            T_total = sum(obj.time_alloc);
        end
    end
end
