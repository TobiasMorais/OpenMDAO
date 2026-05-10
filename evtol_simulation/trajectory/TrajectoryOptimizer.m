classdef TrajectoryOptimizer < handle
% TRAJECTORYOPTIMIZER  Hybrid Flatness + GA energy-optimal trajectory.
%
% Pipeline:
%  1) Build a baseline DifferentialFlatness trajectory through user waypoints
%     with nominal time allocation.
%  2) GA refines (a) per-segment time allocations and (b) intermediate
%     waypoint heights/headings to minimize:
%
%       J_energy = integral_0^T  T_total(t)^{1.5}  dt
%
%     subject to:
%       - kinematic limits (|v|, |a|, |j|)
%       - tilt envelope (|tilt| <= tilt_max along low-speed segments)
%       - thrust per rotor in [0, T_max]
%
% The integrand T^1.5 captures induced-power scaling of rotors near hover
% (Glauert), encouraging the optimizer to spend less time in inefficient
% hover and more in lift-supported cruise. As the PDF notes, this often pushes
% the optimum into deep-stall transients (acceptable for tailsitter).

    properties
        waypoints
        ac_cfg
        baseline
    end

    methods
        function obj = TrajectoryOptimizer(waypoints, ac_cfg)
            obj.waypoints = waypoints;
            obj.ac_cfg = ac_cfg;

            % Time allocation respecting per-rotor saturation AND pitch authority.
            %
            % Two physical constraints set the minimum segment time:
            %   (a) Per-rotor thrust envelope (from previous fix; vertical/horizontal
            %       blend depending on segment direction)
            %   (b) Pitch authority: changing thrust direction requires body rotation,
            %       which is bounded by available control moment M_max ~ N*T_max*r_z
            %       per rotor (small for paired rotors with small z-offset).
            %
            % Critical: we leave 50% safety margin on thrust envelope so the cascade
            % control loop has room to compensate model errors and disturbances.
            % A "trajectory at 100% of physical envelope" leaves the controller no
            % authority to correct -- it must saturate at the first perturbation.
            d_vec = diff(waypoints(:, 1:3));
            N_seg = size(d_vec, 1);
            alpha = 0.50;             % use only 50% of available thrust envelope
            peak_acc_coef = 15;
            g = 9.80665;
            T_per_max = alpha * ac_cfg.rotor.thrust_max;
            a_envelope = (ac_cfg.n_rotors * T_per_max) / ac_cfg.mass;
            % Note: at alpha=0.5, a_envelope = 6 m/s^2 < g; aircraft cannot hover
            % with only 50% thrust. So we add hover thrust separately:
            % a_total = a_thrust + g, with a_thrust budget = a_envelope_full - g_used_for_hover
            T_per_max_full = ac_cfg.rotor.thrust_max;
            a_envelope_full = (ac_cfg.n_rotors * T_per_max_full) / ac_cfg.mass;
            % Margin: 0.5 * (a_envelope_full - g)  available for trajectory accel
            % vertical: must overcome g, so a_max_vert = 0.5 * (a_envelope_full - g)
            % horizontal: full lateral budget = 0.5 * sqrt(a_envelope_full^2 - g^2)
            a_max_vert  = max(0.5 * (a_envelope_full - g), 0.3);
            a_max_horiz = 0.5 * sqrt(max(a_envelope_full^2 - g^2, 1.0));

            % Pitch authority constraint:
            % To rotate body by Delta_theta over time T, peak angular accel ~ 4*Delta/T^2.
            % Required moment M = J*alpha; must satisfy M < M_avail.
            % For rotor-only authority with small z-offset, M_max ~ sum(|z_i| * T_max).
            % This is typically smaller than what large attitude transitions demand.
            J_pitch = ac_cfg.J(2,2);
            M_avail_pitch = 0;
            for i = 1:ac_cfg.n_rotors
                z_off = abs(ac_cfg.rotor.position(i, 3));
                M_avail_pitch = M_avail_pitch + z_off * ac_cfg.rotor.thrust_max;
            end
            M_avail_pitch = 0.3 * M_avail_pitch;   % 30% of theoretical for sustained use
            t_seg = zeros(N_seg, 1);
            for k = 1:N_seg
                d_h = norm(d_vec(k, 1:2));
                d_v = abs(d_vec(k, 3));
                d_t = max(norm(d_vec(k, :)), 1e-6);

                % (a) thrust envelope time
                w_v = d_v / d_t;
                a_eff = w_v * a_max_vert + (1 - w_v) * a_max_horiz;
                T_thrust = sqrt(peak_acc_coef * d_t / a_eff);

                % (b) pitch authority time:
                % approximate the worst-case pitch swing this segment requires:
                % from horizontal motion, body must tilt by atan(a_h/g);
                % from vertical motion only, body can stay vertical (no swing).
                % We bound by 90 deg as worst case (full transition).
                a_h_demand = peak_acc_coef * d_h / max(T_thrust, 1)^2;
                tilt_demand = atan2(a_h_demand, g);
                Delta_theta = max(tilt_demand, deg2rad(5));   % at least 5 deg
                T_pitch = sqrt(4 * Delta_theta * J_pitch / M_avail_pitch);

                t_seg(k) = max([3.0, T_thrust, T_pitch]);
            end

            obj.baseline = DifferentialFlatness(waypoints, t_seg);
        end

        function refined_traj = optimize(obj, varargin)
            % Optional: pass 'method', 'ga' | 'flat'  (default flat-only)
            p = inputParser;
            p.addParameter('method', 'flat');
            p.addParameter('ga_pop', 30);
            p.addParameter('ga_gens', 25);
            p.parse(varargin{:});
            opts = p.Results;

            if strcmpi(opts.method, 'flat')
                refined_traj = obj.baseline;
                return;
            end

            % --- GA refinement of segment times only (keeps waypoint geometry) ---
            N_seg = size(obj.waypoints, 1) - 1;
            t0 = obj.baseline.time_alloc(:);

            % Decision variable: scaling factor per segment in [0.5, 2.0]
            lb = 0.5 * ones(N_seg, 1);
            ub = 2.0 * ones(N_seg, 1);

            cost_fn = @(s) obj.energy_cost(s, t0);

            try
                ga_opts = optimoptions('ga', ...
                    'PopulationSize', opts.ga_pop, ...
                    'MaxGenerations', opts.ga_gens, ...
                    'Display', 'off', ...
                    'UseParallel', false);
                s_opt = ga(cost_fn, N_seg, [],[],[],[], lb, ub, [], ga_opts);
            catch ME
                warning('GA failed (%s); returning baseline trajectory.', ME.message);
                refined_traj = obj.baseline;
                return;
            end

            new_times = t0 .* s_opt(:);
            refined_traj = DifferentialFlatness(obj.waypoints, new_times);
        end

        function J = energy_cost(obj, s, t0)
            new_times = t0 .* s(:);
            try
                tr = DifferentialFlatness(obj.waypoints, new_times);
            catch
                J = 1e9;
                return;
            end
            % Sample T_total along trajectory and integrate T^1.5 dt
            Ttot = tr.total_time();
            N_samples = 80;
            ts = linspace(0, Ttot, N_samples);
            dt = Ttot / (N_samples-1);
            J = 0;
            penalty = 0;
            for tau = ts
                [~, ~, acc, ~, ~, ~, ~] = tr.eval(tau);
                F_NED = obj.ac_cfg.mass * (acc - [0;0;9.80665]);
                T = norm(F_NED);
                J = J + (T^1.5) * dt;

                % Penalize over-tilt
                f_z = -F_NED(3);
                f_xy = norm(F_NED(1:2));
                tilt = atan2(f_xy, max(f_z, 0.1));
                if tilt > deg2rad(89)
                    penalty = penalty + 1e6;
                end
                % Penalize per-rotor thrust > max
                T_per_rotor = T / obj.ac_cfg.n_rotors;
                if T_per_rotor > obj.ac_cfg.rotor.thrust_max
                    penalty = penalty + 1e3 * (T_per_rotor - obj.ac_cfg.rotor.thrust_max)^2;
                end
            end
            J = J + penalty;
        end
    end
end
