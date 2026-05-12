"""Main entry point — runs realistic mission with plots and validation.

Usage:
    python scripts/main.py
"""
from __future__ import annotations

import sys
import os

# Add parent dir to path so 'evtol' is importable when running this script directly
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from evtol.config import aircraft_config, controller_config, simulation_config
from evtol.simulation import run_simulation, build_realistic_mission
from evtol.visualization import plot_telemetry


def main():
    print("=" * 60)
    print("  eVTOL Tailsitter Simulation Framework (Python)")
    print("=" * 60)

    ac = aircraft_config()
    ctrl = controller_config()
    sim = simulation_config()

    print(f"\nMass: {ac['mass']:.0f} kg | Rotors: {ac['n_rotors']} | "
          f"Inner: {ctrl['inner']['type']} | Outer: {ctrl['outer']['type']}")

    # Build realistic mission
    traj, mission_info = build_realistic_mission(ac)
    sim["t_final"] = traj.total_time() + 10
    sim["init_mode"] = mission_info["init_mode"]

    # Run
    log = run_simulation(ac, ctrl, sim, traj)

    # Visualize
    plot_telemetry(log, ac)

    # Final report
    err_norm = np.linalg.norm(log["tracking_err"], axis=0)
    print("\n" + "=" * 60)
    print("  Mission Summary")
    print("=" * 60)
    print(f"  Final position error : {err_norm[-1]:.2f} m")
    print(f"  Max position error   : {np.max(err_norm):.2f} m")
    print(f"  Mean position error  : {np.mean(err_norm):.2f} m")
    print(f"  Final altitude       : {-log['pos_NED'][2, -1]:.1f} m")
    print(f"  Final position N/E   : {log['pos_NED'][0, -1]:.1f} / {log['pos_NED'][1, -1]:.1f} m")
    print("=" * 60)


if __name__ == "__main__":
    main()
