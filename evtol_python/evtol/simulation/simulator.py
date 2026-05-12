"""Main closed-loop simulator.

Cascade architecture per cycle:
    1. Sample reference trajectory at t -> (p_ref, v_ref, a_ref)
    2. Outer PD+FF: f_cmd_NED = K_p*(p_ref-p) + K_v*(v_ref-v) + a_ref - g_NED
    3. force_to_attitude: (qd, T_total, Wd, Wd_dot)
    4. Inner SO(3): virtual moment M_cmd
    5. Allocator: u = [T_1, T_2, T_3, T_4, delta_e, delta_a, delta_r]
    6. Propulsion.step_actuator(T_cmd) + gather forces/moments
    7. Aerodynamics.compute (with slipstream)
    8. Aircraft.dynamics
    9. RK4 integration
"""
from __future__ import annotations

import numpy as np

from evtol.core.aircraft import Aircraft
from evtol.core.quaternion import Quaternion
from evtol.aerodynamics.aerodynamics import Aerodynamics
from evtol.aerodynamics.propulsion import Propulsion
from evtol.control.allocator import ControlAllocator
from evtol.control.so3_controller import AttitudeControllerSO3
from evtol.control.pdff import PositionControllerPDFF
from evtol.control.force_to_attitude import force_to_attitude
from evtol.environment.atmosphere import atmosphere_isa
from evtol.environment.dryden import DrydenWind
from evtol.simulation.rk4 import rk4_step


def _init_log(N: int) -> dict:
    return {
        "t":             np.zeros(N),
        "pos_NED":       np.zeros((3, N)),
        "vel_NED":       np.zeros((3, N)),
        "quat":          np.zeros((4, N)),
        "omega_body":    np.zeros((3, N)),
        "thrust_cmd":    np.zeros((4, N)),
        "thrust_actual": np.zeros((4, N)),
        "elev":          np.zeros(N),
        "ail":           np.zeros(N),
        "rud":           np.zeros(N),
        "pitch_deg":     np.zeros(N),
        "tracking_err":  np.zeros((3, N)),
        "thrust_total":  np.zeros(N),
    }


def _trim_log(log: dict, k: int) -> dict:
    out = {}
    for key, val in log.items():
        if val.ndim == 1:
            out[key] = val[:k]
        else:
            out[key] = val[:, :k]
    return out


def run_simulation(ac_cfg: dict, ctrl_cfg: dict, sim_cfg: dict, traj) -> dict:
    """Run main fixed-step closed-loop simulation.

    Returns log dict with time series of state and controls.
    """
    aircraft = Aircraft(ac_cfg)
    aero = Aerodynamics(ac_cfg)
    prop = Propulsion(ac_cfg, aircraft.rotor_axes)
    alloc = ControlAllocator(ctrl_cfg, ac_cfg, aircraft.rotor_axes)

    # Inner controller (default SO3)
    inner = AttitudeControllerSO3(ctrl_cfg)

    # Outer controller (default PDFF)
    outer_type = ctrl_cfg.get("outer", {}).get("type", "PDFF").upper()
    if outer_type == "PDFF":
        outer = PositionControllerPDFF(ctrl_cfg, ac_cfg["mass"])
    else:
        raise NotImplementedError(f"Outer controller {outer_type} not yet ported")

    wind_model = DrydenWind(sim_cfg["wind"])
    surf_defs = aero.build_surface_strips()

    # Initial state
    init_mode = sim_cfg.get("init_mode", "hover")
    x = aircraft.initial_state(init_mode)

    # Logging buffers
    dt = sim_cfg["dt"]
    N_steps = int(np.floor(sim_cfg["t_final"] / dt)) + 1
    log = _init_log(N_steps)

    # Outer loop trigger
    outer_period = int(round(1 / (ctrl_cfg["f_outer"] * dt)))
    F_cmd_NED = ac_cfg["mass"] * np.array([0.0, 0.0, -9.80665])
    psi_des = 0.0
    prev_attitude_state: dict | None = None

    horizon_pts = ctrl_cfg["nmpc"]["N"] + 1
    horizon_dt = ctrl_cfg["nmpc"]["dt"]

    print(f"Running simulation: {N_steps} steps, t_final={sim_cfg['t_final']:.1f}s")

    for k in range(N_steps):
        t = k * dt

        p_NED = x[0:3]
        v_NED = x[3:6]
        q = x[6:10]
        w_B = x[10:13]
        R_BW = Quaternion.to_rotation_matrix(q)

        # Reference at current time
        p_ref, v_ref, a_ref, _, _, psi_d, _ = traj.eval(min(t, traj.total_time()))
        psi_des = psi_d

        # --- Outer loop ---
        if k % outer_period == 0:
            p_horizon = np.zeros((3, horizon_pts))
            v_horizon = np.zeros((3, horizon_pts))
            for j in range(horizon_pts):
                tj = min(t + j * horizon_dt, traj.total_time())
                pj, vj, _, _, _, _, _ = traj.eval(tj)
                p_horizon[:, j] = pj
                v_horizon[:, j] = vj
            F_specific = outer.compute(p_NED, v_NED, p_horizon, v_horizon, a_ref)
            F_cmd_NED = ac_cfg["mass"] * F_specific

        # --- Force-to-attitude ---
        qd, T_specific, Wd, Wd_dot = force_to_attitude(
            F_cmd_NED / ac_cfg["mass"], psi_des, dt, prev_attitude_state
        )
        prev_attitude_state = {"qd": qd, "Wd": Wd}

        # --- Inner loop (SO3) ---
        Rd = Quaternion.to_rotation_matrix(qd)
        M_cmd = inner.compute(R_BW, w_B, Rd, np.zeros(3), np.zeros(3), ac_cfg["J"], t)

        # Body-frame virtual force: ONLY magnitude along +x_B
        F_virtual_B = np.array([np.linalg.norm(F_cmd_NED), 0.0, 0.0])

        # --- Allocation ---
        v_wind = wind_model.step(max(np.linalg.norm(v_NED), 0.1), max(1.0, -p_NED[2]), dt)
        V_inf_B = R_BW.T @ (v_NED - v_wind)
        rho = atmosphere_isa(-p_NED[2], sim_cfg["env"])
        q_bar = 0.5 * rho * max(np.linalg.norm(V_inf_B), 1.0) ** 2

        nu = np.concatenate([F_virtual_B, M_cmd])
        u_act, _ = alloc.allocate(nu, q_bar)
        T_cmd = u_act[:4]
        delta_e, delta_a, delta_r = u_act[4], u_act[5], u_act[6]

        # --- Propulsion step ---
        prop.step_actuator(T_cmd, dt)
        F_prop_B, M_prop_B, v_i, Omega_signed = prop.compute(V_inf_B, rho)

        # --- Slipstream into surfaces ---
        V_slip = prop.slipstream_velocities(v_i, V_inf_B)

        # --- Aerodynamics ---
        F_aero_B, M_aero_B = aero.compute_forces_moments(V_inf_B, w_B, surf_defs, V_slip, rho)

        # Add surface control contributions
        M_aero_B[1] += (
            -ac_cfg["htail"]["elev_eff"] * ac_cfg["htail"]["area"] * ac_cfg["htail"]["arm"]
        ) * delta_e * q_bar
        M_aero_B[0] += (
            ac_cfg["ail"]["eff"] * ac_cfg["wing"]["area"] * ac_cfg["wing"]["span"]
        ) * delta_a * q_bar
        M_aero_B[2] += (
            -ac_cfg["vtail"]["rud_eff"] * ac_cfg["vtail"]["area"] * ac_cfg["vtail"]["arm"]
        ) * delta_r * q_bar

        # --- Rotor gyroscopic torque ---
        M_gyro_B = aircraft.compute_rotor_gyro(Omega_signed, w_B)

        # --- Integrate ---
        odefun = lambda tt, xx: aircraft.dynamics(xx, F_aero_B, M_aero_B, F_prop_B, M_prop_B, M_gyro_B)
        x = rk4_step(odefun, t, x, dt)

        # --- Log ---
        log["t"][k] = t
        log["pos_NED"][:, k] = p_NED
        log["vel_NED"][:, k] = v_NED
        log["quat"][:, k] = q
        log["omega_body"][:, k] = w_B
        log["thrust_cmd"][:, k] = T_cmd
        log["thrust_actual"][:, k] = ac_cfg["rotor"]["kT"] * prop.Omega_actual ** 2
        log["elev"][k] = delta_e
        log["ail"][k] = delta_a
        log["rud"][k] = delta_r
        eul = Quaternion.to_euler(q)
        log["pitch_deg"][k] = np.rad2deg(eul[1])
        log["tracking_err"][:, k] = p_NED - p_ref
        log["thrust_total"][k] = np.sum(log["thrust_actual"][:, k])

    print(f"Simulation complete: {k+1} steps over {(k+1)*dt:.2f} s")
    return _trim_log(log, k + 1)
