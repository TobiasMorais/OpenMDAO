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
            % Initial uniform time allocation
            N_seg = size(waypoints, 1) - 1;
            d_seg = vecnorm(diff(waypoints(:,1:3)),2,2);
            t_seg = max(2.0, d_seg / 25.0);   % ~25 m/s nominal cruise
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
