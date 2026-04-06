#!/usr/bin/env python3
"""
Plot Artemis II trajectory from Fortran-computed data.

Reads artemis_trajectory.dat (produced by artemis_trajectory.f90)
and creates the same plot as plot_trajectory.py.

Usage:
    python3 plot_artemis_fortran.py                          # uses current time
    python3 plot_artemis_fortran.py 2026-04-06T12:00:00Z     # specific UTC time
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime, timezone


# ── Mission constants ────────────────────────────────────────────────

LAUNCH_TIME = datetime(2026, 4, 1, 22, 35, 0, tzinfo=timezone.utc)
LAUNCH_TS = int(LAUNCH_TIME.timestamp())

def met(days, hours, minutes, seconds=0):
    return days * 86400 + hours * 3600 + minutes * 60 + seconds

PHASES = [
    {"name": "Launch & Ascent",        "color": "#F97316", "start": 0,                    "end": met(0, 3, 23)},
    {"name": "Earth Orbit & Checkout", "color": "#EAB308", "start": met(0, 3, 23),        "end": met(1, 1, 8, 42)},
    {"name": "Trans-Lunar Injection",  "color": "#EC4899", "start": met(1, 1, 8, 42),     "end": met(1, 1, 38, 42)},
    {"name": "Outbound Coast",         "color": "#4ADE80", "start": met(1, 1, 38, 42),    "end": met(4, 18, 0)},
    {"name": "Lunar Flyby",            "color": "#A855F7", "start": met(4, 18, 0),        "end": met(5, 18, 53)},
    {"name": "Return Coast",           "color": "#60A5FA", "start": met(5, 18, 53),       "end": met(9, 1, 9)},
    {"name": "Re-entry & Splashdown",  "color": "#F97316", "start": met(9, 1, 9),         "end": met(9, 2, 0)},
]


# ── Load Fortran output ─────────────────────────────────────────────

def load_trajectory(path="artemis_trajectory.dat"):
    """Load trajectory data produced by the Fortran program."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            rows.append({
                "timestamp":       int(parts[0]),
                "met_s":           float(parts[1]),
                "x":               float(parts[2]),
                "y":               float(parts[3]),
                "z":               float(parts[4]),
                "moon_x":          float(parts[5]),
                "moon_y":          float(parts[6]),
                "moon_z":          float(parts[7]),
                "earth_distance_km": float(parts[8]),
                "moon_distance_km":  float(parts[9]),
                "velocity_km_s":   float(parts[10]),
                "range_rate_km_s": float(parts[11]),
            })
    print(f"Loaded {len(rows)} trajectory points from {path}")
    return rows


# ── Helpers ──────────────────────────────────────────────────────────

def get_phase_color(timestamp):
    met_s = timestamp - LAUNCH_TS
    for phase in PHASES:
        if phase["start"] <= met_s < phase["end"]:
            return phase["color"]
    return "#888888"


def interpolate(rows, timestamp):
    if not rows:
        return None
    idx = 0
    for i, r in enumerate(rows):
        if r["timestamp"] <= timestamp:
            idx = i
    if idx < len(rows) - 1 and timestamp >= rows[idx]["timestamp"]:
        a, b = rows[idx], rows[idx + 1]
        f = min(1.0, (timestamp - a["timestamp"]) / (b["timestamp"] - a["timestamp"]))
        return {k: a[k] + (b[k] - a[k]) * f for k in a}
    return rows[idx]


# ── Plot ─────────────────────────────────────────────────────────────

def plot_trajectory(rows, target_time):
    fig, ax = plt.subplots(figsize=(12, 10))
    fig.patch.set_facecolor("#0f172a")
    ax.set_facecolor("#0f172a")

    # Trajectory segments colored by mission phase
    for i in range(len(rows) - 1):
        color = get_phase_color(rows[i]["timestamp"])
        ax.plot(
            [rows[i]["x"], rows[i + 1]["x"]],
            [rows[i]["y"], rows[i + 1]["y"]],
            color=color, linewidth=0.8, alpha=0.85,
        )

    # Moon trajectory (dashed)
    ax.plot(
        [r["moon_x"] for r in rows],
        [r["moon_y"] for r in rows],
        color="#555555", linewidth=0.5, alpha=0.5, linestyle="--",
    )

    # Earth at origin
    ax.plot(0, 0, "o", color="#3b82f6", markersize=14, zorder=10)
    ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(12, -5),
                color="white", fontsize=10, fontweight="bold")

    # Current position marker
    target_ts = int(target_time.timestamp())
    interp = interpolate(rows, target_ts)

    if interp:
        # Moon
        ax.plot(interp["moon_x"], interp["moon_y"], "o", color="#d1d5db",
                markersize=10, zorder=10)
        ax.annotate("Moon", (interp["moon_x"], interp["moon_y"]),
                    textcoords="offset points", xytext=(10, -5),
                    color="#d1d5db", fontsize=10, fontweight="bold")

        # Orion
        ax.plot(interp["x"], interp["y"], "o", color="#ffffff",
                markersize=8, zorder=11, markeredgecolor="#facc15",
                markeredgewidth=2)
        met_s = target_ts - LAUNCH_TS
        met_days = met_s // 86400
        met_hours = (met_s % 86400) // 3600
        met_mins = (met_s % 3600) // 60
        label = (
            f"Orion @ {target_time.strftime('%Y-%m-%d %H:%M UTC')}\n"
            f"MET T+{met_days}d {met_hours:02d}h {met_mins:02d}m\n"
            f"Earth: {interp['earth_distance_km']:,.0f} km\n"
            f"Moon: {interp['moon_distance_km']:,.0f} km\n"
            f"Speed: {interp['velocity_km_s']:.2f} km/s"
        )
        ax.annotate(label, (interp["x"], interp["y"]),
                    textcoords="offset points", xytext=(15, 10),
                    color="#facc15", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.4", fc="#1e293b", ec="#facc15", alpha=0.9))

    # Phase legend
    legend_patches = [mpatches.Patch(color=p["color"], label=p["name"]) for p in PHASES]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=8,
              facecolor="#1e293b", edgecolor="#334155", labelcolor="white")

    ax.set_xlabel("X (km, J2000 Earth-centered)", color="white", fontsize=11)
    ax.set_ylabel("Y (km, J2000 Earth-centered)", color="white", fontsize=11)
    ax.set_title("Artemis II — Orion Trajectory (Fortran-computed)", color="white",
                 fontsize=15, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("#334155")
    ax.set_aspect("equal")
    ax.grid(True, color="#1e293b", linewidth=0.5)

    out_path = "artemis_trajectory_fortran.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close()
    print(f"\nSaved plot to {out_path}")


# ── Entry point ──────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) > 1:
        target_time = datetime.fromisoformat(sys.argv[1].replace("Z", "+00:00"))
    else:
        target_time = datetime.now(timezone.utc)

    print(f"Target time: {target_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")

    rows = load_trajectory()
    plot_trajectory(rows, target_time)
