"""Physical sanity checks for the eVTOL model.

Validates 20+ first-principles indicators against typical ranges for
this class of aircraft (Joby S4, Beta Alia).
"""
from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from evtol.config import aircraft_config, controller_config
from evtol.core import Aircraft


def _print_check(name: str, value: str, ok: bool | None, note: str = "") -> None:
    tag = "[OK]  " if ok is True else ("[WARN]" if ok is False else "[INFO]")
    print(f"  {tag} {name:<32} = {value:<30} | {note}")


def main():
    ac = aircraft_config()
    ctrl = controller_config()
    g = 9.80665
    rho = 1.225
    W = ac["mass"] * g

    print("\n" + "=" * 67)
    print("  eVTOL Tailsitter — Physical Sanity Checks")
    print("=" * 67 + "\n")

    _print_check("Mass (m)", f"{ac['mass']:.0f} kg", None)
    _print_check("Weight (W = m*g)", f"{W:.0f} N", None)

    # T/W
    T_avail = ac["n_rotors"] * ac["rotor"]["thrust_max"]
    TW = T_avail / W
    _print_check("T/W ratio", f"{TW:.2f}", 1.10 < TW < 2.0, "Tailsitter typical 1.2-1.6")

    # Hover throttle
    T_hover = W / ac["n_rotors"]
    throttle = T_hover / ac["rotor"]["thrust_max"]
    _print_check("Hover throttle", f"{100*throttle:.1f}%", 0.5 < throttle < 0.95,
                 "Healthy: 60-85%")

    # Wing loading
    WL_kg = (W / ac["wing"]["area"]) / g
    _print_check("Wing loading W/S", f"{WL_kg:.1f} kg/m^2",
                 50 < WL_kg < 250, "Light eVTOL 80-180")

    # Disk loading
    A_total = ac["n_rotors"] * ac["rotor"]["disk_area"]
    DL_kg = (W / A_total) / g
    _print_check("Disk loading T/A", f"{DL_kg:.1f} kg/m^2",
                 50 < DL_kg < 800, "Open rotor 100-500")

    # Hover induced velocity
    v_i = np.sqrt(T_hover / (2 * rho * ac["rotor"]["disk_area"]))
    _print_check("Hover induced velocity", f"{v_i:.2f} m/s", 5 < v_i < 60,
                 "Open rotor 10-30 m/s")

    # Hover power (with FM=0.7)
    P_total_kW = ac["n_rotors"] * T_hover * v_i / 0.7 / 1000
    _print_check("Hover power total", f"{P_total_kW:.0f} kW",
                 50 < P_total_kW < 1500, "Light eVTOL 100-500 kW")
    PW = ac["mass"] / P_total_kW
    _print_check("Power loading (kg/kW)", f"{PW:.2f} kg/kW",
                 1.5 < PW < 8.0, "Healthy 2-5")

    # Stall speed
    V_stall = np.sqrt(2 * W / (rho * ac["wing"]["area"] * ac["wing"]["CL_max"]))
    _print_check("Stall speed V_s", f"{V_stall:.1f} m/s",
                 15 < V_stall < 50, "Tailsitter 25-45")

    # Cruise at 1.3*V_s
    V_cruise = 1.3 * V_stall
    CL_cruise = 2 * W / (rho * ac["wing"]["area"] * V_cruise ** 2)
    alpha_cruise = CL_cruise / ac["wing"]["CL_alpha"] + ac["wing"]["alpha0"]
    K = 1 / (np.pi * ac["wing"]["AR"] * ac["wing"]["e"])
    CD_cruise = ac["wing"]["CD0"] + K * CL_cruise ** 2
    LD = CL_cruise / CD_cruise
    _print_check("Cruise speed (1.3 V_s)", f"{V_cruise:.1f} m/s",
                 25 < V_cruise < 60, "Reasonable for tailsitter")
    _print_check("Cruise CL", f"{CL_cruise:.3f}", 0.3 < CL_cruise < ac["wing"]["CL_max"],
                 f"Below CL_max={ac['wing']['CL_max']:.2f}")
    _print_check("Cruise alpha", f"{np.rad2deg(alpha_cruise):.2f} deg",
                 abs(np.rad2deg(alpha_cruise)) < np.rad2deg(ac["wing"]["alpha_stall"]),
                 "Below stall AoA")
    _print_check("Cruise L/D", f"{LD:.1f}", 6 < LD < 25, "Light fixed-wing 8-15")

    # Max L/D point
    CL_LDmax = np.sqrt(ac["wing"]["CD0"] / K)
    V_LDmax = np.sqrt(2 * W / (rho * ac["wing"]["area"] * CL_LDmax))
    LD_max = 0.5 / np.sqrt(ac["wing"]["CD0"] * K)
    _print_check("V_LDmax", f"{V_LDmax:.1f} m/s",
                 25 < V_LDmax < 60, f"alpha_LDmax = {np.rad2deg(CL_LDmax/ac['wing']['CL_alpha']):.1f} deg")
    _print_check("L/D_max", f"{LD_max:.1f}", 8 < LD_max < 25, "")

    # Bandwidth separation
    tau_motor = ac["rotor"]["tau_motor"]
    omega_n_inner = np.sqrt(ctrl["so3"]["kR"][1, 1] / ac["J"][1, 1])
    tau_inner = 1 / (omega_n_inner * 0.7)
    # PDFF outer: omega_n = 1 rad/s
    omega_n_outer = 1.0
    tau_outer = 1 / (omega_n_outer * 0.85)
    r_im = tau_inner / tau_motor
    r_oi = tau_outer / tau_inner
    _print_check("Motor tau", f"{tau_motor*1000:.0f} ms",
                 tau_motor > 0.01 and tau_motor < 0.2, "Typical 30-100 ms")
    _print_check("Inner SO3 tau", f"{tau_inner*1000:.0f} ms",
                 tau_inner > tau_motor * 1.5,
                 f"Separation inner/motor = {r_im:.1f}x (need >1.5)")
    _print_check("Outer PDFF tau", f"{tau_outer*1000:.0f} ms",
                 tau_outer > tau_inner * 1.5,
                 f"Separation outer/inner = {r_oi:.1f}x (need >1.5)")

    # Inertia ratios
    J_xx = ac["J"][0, 0]; J_yy = ac["J"][1, 1]; J_zz = ac["J"][2, 2]; J_xz = abs(ac["J"][0, 2])
    _print_check("J_yy/J_xx (pitch vs roll)", f"{J_yy/J_xx:.2f}",
                 1.0 < J_yy/J_xx < 4.0, "Tailsitter has more pitch inertia")
    _print_check("J_xz/J_yy (coupling)", f"{J_xz/J_yy:.3f}", J_xz/J_yy < 0.20,
                 "Should be small")

    # Free fall and hover trim
    aircraft = Aircraft(ac)
    x0 = np.concatenate([np.zeros(6), np.array([1, 0, 0, 0]), np.zeros(3)])
    xdot = aircraft.dynamics(x0, np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3))
    _print_check("Free-fall acceleration", f"|a|={np.linalg.norm(xdot[3:6]):.4f}",
                 abs(np.linalg.norm(xdot[3:6]) - g) < 1e-9, "No forces -> a = g")

    x_h = aircraft.initial_state("hover")
    F_prop_B = np.array([W, 0, 0])
    xdot_h = aircraft.dynamics(x_h, np.zeros(3), np.zeros(3), F_prop_B, np.zeros(3), np.zeros(3))
    _print_check("Hover trim (T = m*g)", f"|a|={np.linalg.norm(xdot_h[3:6]):.3e}",
                 np.linalg.norm(xdot_h[3:6]) < 1e-6, "Should give zero net accel")

    print("\n" + "=" * 67)
    print("  Summary: physical sanity check complete")
    print("=" * 67 + "\n")


if __name__ == "__main__":
    main()
