"""Aerospace-grade validation plots from simulation log."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_telemetry(log: dict, ac_cfg: dict) -> None:
    """Generate 8 standard plots: 3D trajectory, NED pos/vel, rates, pitch, motor, error."""
    t = log["t"]

    # 1. 3D Trajectory
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(log["pos_NED"][1], log["pos_NED"][0], -log["pos_NED"][2], "b-", linewidth=1.5)
    ax.plot(log["pos_NED"][1, 0], log["pos_NED"][0, 0], -log["pos_NED"][2, 0],
            "go", markersize=10, label="start")
    ax.plot(log["pos_NED"][1, -1], log["pos_NED"][0, -1], -log["pos_NED"][2, -1],
            "rs", markersize=10, label="end")
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.set_zlabel("Altitude [m]")
    ax.set_title("3D Trajectory (start=green, end=red)")
    ax.legend()

    # 2. Position NED
    fig, axes = plt.subplots(3, 1, figsize=(10, 8))
    axes[0].plot(t, log["pos_NED"][0], "b")
    axes[0].set_ylabel("North [m]"); axes[0].grid()
    axes[0].set_title("Position (NED frame)")
    axes[1].plot(t, log["pos_NED"][1], "b")
    axes[1].set_ylabel("East [m]"); axes[1].grid()
    axes[2].plot(t, -log["pos_NED"][2], "b")
    axes[2].set_ylabel("Altitude [m]"); axes[2].set_xlabel("t [s]"); axes[2].grid()

    # 3. Velocity NED
    fig, axes = plt.subplots(3, 1, figsize=(10, 8))
    axes[0].plot(t, log["vel_NED"][0], "b"); axes[0].set_ylabel("V_N [m/s]"); axes[0].grid()
    axes[1].plot(t, log["vel_NED"][1], "b"); axes[1].set_ylabel("V_E [m/s]"); axes[1].grid()
    axes[2].plot(t, log["vel_NED"][2], "b"); axes[2].set_ylabel("V_D [m/s]"); axes[2].set_xlabel("t [s]"); axes[2].grid()

    # 4. Body rates
    plt.figure(figsize=(10, 6))
    plt.plot(t, np.rad2deg(log["omega_body"][0]), "r", label="p")
    plt.plot(t, np.rad2deg(log["omega_body"][1]), "g", label="q")
    plt.plot(t, np.rad2deg(log["omega_body"][2]), "b", label="r")
    plt.grid(); plt.legend(); plt.ylabel("Body rates [deg/s]"); plt.xlabel("t [s]")
    plt.title("Angular rates in body frame")

    # 5. Pitch transition
    plt.figure(figsize=(10, 6))
    plt.plot(t, log["pitch_deg"], "b-", linewidth=1.5)
    plt.axhline(90, color="k", linestyle="--", label="Hover (90°)")
    plt.axhline(0,  color="k", linestyle="--", label="Cruise (0°)")
    plt.grid(); plt.legend(); plt.xlabel("t [s]"); plt.ylabel(r"$\theta$ [deg]")
    plt.title("Pitch angle vs time")

    # 6. Per-rotor thrust
    plt.figure(figsize=(10, 6))
    plt.plot(t, log["thrust_actual"][0], "r", label="Rotor 1 (R, upper)")
    plt.plot(t, log["thrust_actual"][1], "g", label="Rotor 2 (R, lower)")
    plt.plot(t, log["thrust_actual"][2], "b", label="Rotor 3 (L, upper)")
    plt.plot(t, log["thrust_actual"][3], "m", label="Rotor 4 (L, lower)")
    plt.grid(); plt.legend(); plt.ylabel("Thrust [N]"); plt.xlabel("t [s]")
    plt.title("Per-rotor actual thrust")

    # 7. Differential motor effort
    plt.figure(figsize=(10, 6))
    avg = np.mean(log["thrust_actual"], axis=0)
    diff = log["thrust_actual"] - avg
    for i, label in enumerate(["R1", "R2", "R3", "R4"]):
        plt.plot(t, diff[i], linewidth=1.0, label=label)
    plt.grid(); plt.legend(); plt.ylabel(r"$\Delta T_i$ [N]"); plt.xlabel("t [s]")
    plt.title(r"Differential thrust per rotor ($T_i - \bar T$)")

    # 8. Tracking error
    plt.figure(figsize=(10, 6))
    err = np.linalg.norm(log["tracking_err"], axis=0)
    plt.plot(t, err, "b-", linewidth=1.5)
    plt.grid(); plt.ylabel(r"$\|p - p_{ref}\|$ [m]"); plt.xlabel("t [s]")
    plt.title("Position tracking error")

    print(f"\nMaximum tracking error: {np.max(err):.2f} m")
    energy = np.trapz(log["thrust_total"] ** 1.5, t)
    print(f"Total energy proxy J = int T^1.5 dt = {energy:.3e}")

    plt.show()
