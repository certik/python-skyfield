#!/usr/bin/env python3
"""
Artemis II Trajectory Plotter

Downloads spacecraft (Orion) and Moon ephemeris from JPL Horizons,
then plots the Earth-centered 2D trajectory with the current position.

Usage:
    python3 plot_trajectory.py                          # uses current time
    python3 plot_trajectory.py 2026-04-06T12:00:00Z     # specific UTC time
"""

import sys
import math
import json
import os
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime, timezone

# ── Mission constants ────────────────────────────────────────────────

ORION_ID = "-1024"
MOON_ID = "301"
EARTH_CENTER = "500@399"

LAUNCH_TIME = datetime(2026, 4, 1, 22, 35, 0, tzinfo=timezone.utc)
EPHEMERIS_START = datetime(2026, 4, 2, 2, 0, 0, tzinfo=timezone.utc)
MISSION_END = datetime(2026, 4, 10, 23, 0, 0, tzinfo=timezone.utc)
STEP_SIZE = "5 m"  # 5-minute intervals
DATA_FILE = "artemis_ephemeris.json"


# ── Mission phases (MET-based, same as milestones.ts) ────────────────

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

LAUNCH_TS = int(LAUNCH_TIME.timestamp())


# ── JPL Horizons API ─────────────────────────────────────────────────

def fetch_vectors(command, center, start_time, stop_time, step_size):
    """Fetch vector ephemeris from JPL Horizons."""
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": f"'{center}'",
        "START_TIME": f"'{start_time}'",
        "STOP_TIME": f"'{stop_time}'",
        "STEP_SIZE": f"'{step_size}'",
        "OUT_UNITS": "'KM-S'",
        "REF_SYSTEM": "'J2000'",
        "VEC_TABLE": "'3'",
        "CSV_FORMAT": "'YES'",
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    print(f"  Fetching {command} from Horizons...")
    resp = requests.get(url, params=params, timeout=120)
    resp.raise_for_status()
    return resp.json()


def parse_records(result):
    """Parse full records (position + velocity + range) from Horizons response."""
    text = result["result"]
    soe = text.index("$$SOE")
    eoe = text.index("$$EOE")
    records = []
    for line in text[soe + 5:eoe].strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        p = [s.strip() for s in line.split(",")]
        jd = float(p[0])
        ts = int((jd - 2440587.5) * 86400)
        records.append({
            "timestamp": ts,
            "x": float(p[2]), "y": float(p[3]), "z": float(p[4]),
            "vx": float(p[5]), "vy": float(p[6]), "vz": float(p[7]),
            "range": float(p[9]), "range_rate": float(p[10]),
        })
    return records


def parse_pos_records(result):
    """Parse position-only records from Horizons response."""
    text = result["result"]
    soe = text.index("$$SOE")
    eoe = text.index("$$EOE")
    records = []
    for line in text[soe + 5:eoe].strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        p = [s.strip() for s in line.split(",")]
        jd = float(p[0])
        ts = int((jd - 2440587.5) * 86400)
        records.append({"timestamp": ts, "x": float(p[2]), "y": float(p[3]), "z": float(p[4])})
    return records


# ── Main logic ───────────────────────────────────────────────────────

def download_ephemeris():
    """Download Orion and Moon ephemeris for the full mission."""
    start_str = EPHEMERIS_START.strftime("%Y-%m-%d %H:%M:%S")
    stop_str = MISSION_END.strftime("%Y-%m-%d %H:%M:%S")

    orion_raw = fetch_vectors(ORION_ID, EARTH_CENTER, start_str, stop_str, STEP_SIZE)
    moon_raw = fetch_vectors(MOON_ID, EARTH_CENTER, start_str, stop_str, STEP_SIZE)

    orion = parse_records(orion_raw)
    moon = parse_pos_records(moon_raw)
    moon_by_ts = {m["timestamp"]: m for m in moon}

    rows = []
    for o in orion:
        m = moon_by_ts.get(o["timestamp"])
        if not m:
            continue
        dx, dy, dz = o["x"] - m["x"], o["y"] - m["y"], o["z"] - m["z"]
        rows.append({
            "timestamp": o["timestamp"],
            "earth_distance_km": o["range"],
            "moon_distance_km": math.sqrt(dx * dx + dy * dy + dz * dz),
            "velocity_km_s": math.sqrt(o["vx"] ** 2 + o["vy"] ** 2 + o["vz"] ** 2),
            "range_rate_km_s": o["range_rate"],
            "x": o["x"], "y": o["y"], "z": o["z"],
            "moon_x": m["x"], "moon_y": m["y"], "moon_z": m["z"],
        })
    return rows


def get_phase_color(timestamp):
    """Return the phase color for a given unix timestamp."""
    met_s = timestamp - LAUNCH_TS
    for phase in PHASES:
        if phase["start"] <= met_s < phase["end"]:
            return phase["color"]
    return "#888888"


def interpolate(rows, timestamp):
    """Linear interpolation of a telemetry row at exact timestamp."""
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


def plot_trajectory(rows, target_time):
    """Plot 2D trajectory (X-Y plane, Earth-centered) with phase colors."""
    fig, ax = plt.subplots(figsize=(12, 10))
    fig.patch.set_facecolor("#0f172a")
    ax.set_facecolor("#0f172a")

    # Draw trajectory segments colored by phase
    for i in range(len(rows) - 1):
        color = get_phase_color(rows[i]["timestamp"])
        ax.plot(
            [rows[i]["x"], rows[i + 1]["x"]],
            [rows[i]["y"], rows[i + 1]["y"]],
            color=color, linewidth=0.8, alpha=0.85,
        )

    # Draw Moon trajectory
    ax.plot(
        [r["moon_x"] for r in rows],
        [r["moon_y"] for r in rows],
        color="#555555", linewidth=0.5, alpha=0.5, linestyle="--",
    )

    # Earth at origin
    ax.plot(0, 0, "o", color="#3b82f6", markersize=14, zorder=10)
    ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(12, -5),
                color="white", fontsize=10, fontweight="bold")

    # Moon position at target time
    target_ts = int(target_time.timestamp())
    interp = interpolate(rows, target_ts)

    if interp:
        ax.plot(interp["moon_x"], interp["moon_y"], "o", color="#d1d5db",
                markersize=10, zorder=10)
        ax.annotate("Moon", (interp["moon_x"], interp["moon_y"]),
                    textcoords="offset points", xytext=(10, -5),
                    color="#d1d5db", fontsize=10, fontweight="bold")

        # Orion current position
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
    ax.set_title("Artemis II — Orion Trajectory", color="white", fontsize=15, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("#334155")
    ax.set_aspect("equal")
    ax.grid(True, color="#1e293b", linewidth=0.5)

    out_path = "artemis_trajectory.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close()
    print(f"\nSaved plot to {out_path}")


def save_ephemeris(rows, path=DATA_FILE):
    """Save ephemeris data to JSON."""
    with open(path, "w") as f:
        json.dump(rows, f)
    print(f"Saved {len(rows)} data points to {path}")


def load_ephemeris(path=DATA_FILE):
    """Load cached ephemeris data from JSON, or return None."""
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = json.load(f)
    print(f"Loaded {len(rows)} data points from {path}")
    return rows


# ── Entry point ──────────────────────────────────────────────────────

if __name__ == "__main__":
    # Parse optional target time
    if len(sys.argv) > 1:
        target_time = datetime.fromisoformat(sys.argv[1].replace("Z", "+00:00"))
    else:
        target_time = datetime.now(timezone.utc)

    print(f"Target time: {target_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")

    rows = load_ephemeris()
    if rows is None:
        print("Downloading ephemeris from JPL Horizons...")
        rows = download_ephemeris()
        print(f"  Got {len(rows)} data points")
        save_ephemeris(rows)

    plot_trajectory(rows, target_time)
