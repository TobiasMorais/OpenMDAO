function log = run_simulation(ac_cfg, ctrl_cfg, sim_cfg, traj)
% RUN_SIMULATION  Main fixed-step closed-loop simulation.
%
% Cascade architecture per cycle:
%   1) Sample reference trajectory at t -> (p_ref, v_ref, F_des_ff)
%   2) Outer NMPC: F_cmd_NED = NMPC( state, ref_horizon )
%   3) force_to_attitude:  (qd, T_total, Wd, Wd_dot)
%   4) Inner SO(3) or INDI: virtual moment M_cmd
%   5) Allocator: u = [T_1..T_4, delta_e, delta_a, delta_r]
%   6) Propulsion.step_actuator(T_cmd) and gather forces/moments
%   7) Aerodynamics.compute (with slipstream)
%   8) Aircraft.dynamics
%   9) RK4 integration
%
% Loop runs at sim_cfg.dt (typically 500 Hz).
% Outer loop runs every (1/f_outer)/dt steps.

    aircraft = Aircraft(ac_cfg);
    aero     = Aerodynamics(ac_cfg);
    prop     = Propulsion(ac_cfg, aircraft.rotor_axes);
    alloc    = ControlAllocator(ctrl_cfg, ac_cfg, aircraft.rotor_axes);

    % Inner controller selection
    switch upper(ctrl_cfg.inner.type)
        case 'SO3'
            inner = AttitudeControllerSO3(ctrl_cfg);
        case 'INDI'
            G_eff = alloc.attitude_effectiveness();
            inner = AttitudeControllerINDI(ctrl_cfg, G_eff, ac_cfg.n_rotors);
        otherwise
            error('Unknown inner controller: %s', ctrl_cfg.inner.type);
    end

    outer = PositionControllerNMPC(ctrl_cfg, ac_cfg.mass);
    wind_model = DrydenWind(sim_cfg.wind);
    surf_defs  = aero.build_surface_strips();

    % Initial state (override via sim_cfg.init_mode if provided)
    init_mode = 'hover';
    if isfield(sim_cfg, 'init_mode') && ~isempty(sim_cfg.init_mode)
        init_mode = sim_cfg.init_mode;
    end
    x = aircraft.initial_state(init_mode);

    % Logging buffers — run for t_final regardless of trajectory length.
    % After traj.total_time() the reference is held at the last waypoint.
    N_steps = floor(sim_cfg.t_final / sim_cfg.dt) + 1;
    log = init_log(N_steps);

    % Outer loop trigger
    outer_period = round(1 / (ctrl_cfg.f_outer * sim_cfg.dt));
    F_cmd_NED = ac_cfg.mass * [0;0;-9.80665];   % init: hold against gravity
    psi_des = 0.0;
    prev_attitude_state = struct('qd', [], 'Wd', []);

    % Outer-step ref preview (precompute only when called)
    horizon_pts = ctrl_cfg.nmpc.N + 1;
    horizon_dt  = ctrl_cfg.nmpc.dt;

    for k = 1:N_steps
        t = (k-1) * sim_cfg.dt;

        p_NED = x(1:3);
        v_NED = x(4:6);
        q     = x(7:10);
        w_B   = x(11:13);

        R_BW = quat_utils('toR', q);

        % --- Reference trajectory sample ---
        [p_ref, v_ref, ~, ~, ~, psi_d, ~] = traj.eval(min(t, traj.total_time()));
        psi_des = psi_d;

        % --- Outer loop (NMPC) ---
        if mod(k-1, outer_period) == 0
            % Build short reference horizon by sampling traj
            p_horizon = zeros(3, horizon_pts);
            v_horizon = zeros(3, horizon_pts);
            for j = 0:horizon_pts-1
                tj = min(t + j*horizon_dt, traj.total_time());
                [pj, vj, ~, ~, ~, ~, ~] = traj.eval(tj);
                p_horizon(:, j+1) = pj;
                v_horizon(:, j+1) = vj;
            end
            F_specific = outer.compute(p_NED, v_NED, p_horizon, v_horizon);
            F_cmd_NED = ac_cfg.mass * F_specific;
        end

        % --- Force-to-attitude mapping ---
        % NMPC outputs specific thrust f_cmd = F_cmd_NED/m (already gravity-decoupled
        % since predictive model is v_dot = f_cmd + g_NED). Thrust direction = f_cmd.
        [qd, T_specific, Wd, Wd_dot] = force_to_attitude( ...
            F_cmd_NED / ac_cfg.mass, psi_des, sim_cfg.dt, prev_attitude_state);
        T_total_cmd = T_specific * ac_cfg.mass;
        prev_attitude_state.qd = qd;
        prev_attitude_state.Wd = Wd;

        % --- Inner attitude loop ---
        % NOTE: We pass ZERO Wd, Wd_dot to inner controllers because the FD-
        % based estimates from force_to_attitude spike at every NMPC update
        % (qd jumps when F_cmd_NED jumps every outer_period steps). Pure-
        % feedback SO(3) is slightly slower but unconditionally stable, while
        % the FD-amplified FF can drive divergence. If smoother qd profiles
        % are available (e.g., from differential flatness recover_states), Wd
        % can be re-enabled with care.
        switch upper(ctrl_cfg.inner.type)
            case 'SO3'
                Rd = quat_utils('toR', qd);
                M_cmd = inner.compute(R_BW, w_B, Rd, zeros(3,1), zeros(3,1), ac_cfg.J, t);
                % Body-frame virtual force: only the magnitude along +x_B
                % (which is what the rotors physically produce). The body is
                % rotated to align x_B with F_cmd_NED by the SO(3) loop;
                % giving the allocator off-axis components would force it to
                % use surfaces/diff-thrust to fight a misalignment that SO(3)
                % is already correcting -- causes loop instability.
                F_virtual_B = [norm(F_cmd_NED); 0; 0];
            case 'INDI'
                u_meas = [prop.Omega_actual.^2 * ac_cfg.rotor.kT; 0; 0; 0];
                [Du, ~] = inner.compute(q, w_B, qd, zeros(3,1), u_meas, t);
                M_cmd = ac_cfg.J * (Du(1:3));   % approximate: pretend Du is angular accel demand
                G_att = alloc.attitude_effectiveness();
                M_cmd = ac_cfg.J * (G_att * Du(1:ac_cfg.n_rotors));
                F_virtual_B = [norm(F_cmd_NED); 0; 0];
        end

        % --- Allocation ---
        nu = [F_virtual_B; M_cmd];
        V_inf_B = R_BW' * (v_NED - wind_step(wind_model, v_NED, p_NED, sim_cfg.dt));
        rho = atmosphere_isa(-p_NED(3), sim_cfg.env);
        q_bar = 0.5 * rho * max(norm(V_inf_B), 1.0)^2;
        [u_act, ~] = alloc.allocate(nu, q_bar);

        T_cmd_per_rotor = u_act(1:4);
        delta_e = u_act(5);  delta_a = u_act(6);  delta_r = u_act(7);

        % --- Propulsion step ---
        prop.step_actuator(T_cmd_per_rotor, sim_cfg.dt);
        [F_prop_B, M_prop_B, v_i, Omega_signed] = prop.compute(V_inf_B, rho);

        % --- Slipstream into surfaces ---
        V_slip = prop.slipstream_velocities(v_i, V_inf_B);

        % --- Aerodynamics ---
        [F_aero_B, M_aero_B] = aero.compute_forces_moments( ...
            V_inf_B, w_B, surf_defs, V_slip, rho);

        % Add control surface direct contribution to moments
        M_aero_B(2) = M_aero_B(2) + (-ac_cfg.htail.elev_eff * ac_cfg.htail.area * ac_cfg.htail.arm) * delta_e * q_bar;
        M_aero_B(1) = M_aero_B(1) + ac_cfg.ail.eff * ac_cfg.wing.area * ac_cfg.wing.span * delta_a * q_bar;
        M_aero_B(3) = M_aero_B(3) + (-ac_cfg.vtail.rud_eff * ac_cfg.vtail.area * ac_cfg.vtail.arm) * delta_r * q_bar;

        % --- Rotor gyroscopic torque ---
        M_gyro_B = aircraft.compute_rotor_gyro(Omega_signed, w_B);

        % --- Pack ode and integrate ---
        odefun = @(tt, xx) aircraft.dynamics(xx, F_aero_B, M_aero_B, F_prop_B, M_prop_B, M_gyro_B);
        x = rk4_step(odefun, t, x, sim_cfg.dt);

        % --- Logging ---
        log.t(k) = t;
        log.pos_NED(:,k)    = p_NED;
        log.vel_NED(:,k)    = v_NED;
        log.quat(:,k)       = q;
        log.omega_body(:,k) = w_B;
        log.thrust_cmd(:,k) = T_cmd_per_rotor;
        log.thrust_actual(:,k) = ac_cfg.rotor.kT * prop.Omega_actual.^2;
        log.elev(k) = delta_e; log.ail(k) = delta_a; log.rud(k) = delta_r;
        eul = quat_utils('toEuler', q);
        log.pitch_deg(k) = rad2deg(eul(2));
        % Tracking error (current pos minus reference)
        log.tracking_err(:,k) = p_NED - p_ref;
        % Energy proxy
        T_total = sum(log.thrust_actual(:,k));
        log.thrust_total(k) = T_total;
    end

    log.t = log.t(1:k);
    log = trim_log(log, k);
    fprintf('Simulation complete: %d steps over %.2f s\n', k, log.t(end));
end

function v_wind_NED = wind_step(wind_model, v_NED, p_NED, dt)
    Vair = norm(v_NED);
    h    = max(1.0, -p_NED(3));
    v_wind_NED = wind_model.step(Vair, h, dt);
end

function log = init_log(N)
    log.t = zeros(1,N);
    log.pos_NED = zeros(3,N);
    log.vel_NED = zeros(3,N);
    log.quat = zeros(4,N);
    log.omega_body = zeros(3,N);
    log.thrust_cmd = zeros(4,N);
    log.thrust_actual = zeros(4,N);
    log.elev = zeros(1,N);
    log.ail = zeros(1,N);
    log.rud = zeros(1,N);
    log.pitch_deg = zeros(1,N);
    log.tracking_err = zeros(3,N);
    log.thrust_total = zeros(1,N);
end

function log = trim_log(log, k)
    f = fieldnames(log);
    for i = 1:numel(f)
        v = log.(f{i});
        if size(v,2) >= k
            log.(f{i}) = v(:, 1:k);
        end
    end
end
