#!/usr/bin/env python3
"""
Convert artemis_ephemeris.json to a simple text format readable by Fortran.

Input:  artemis_ephemeris.json  (from plot_trajectory.py / JPL Horizons)
Output: artemis_orion.dat       (text table for Fortran)

Run:  python3 convert_artemis_data.py
"""

import json
from datetime import datetime, timezone

DATA_FILE = "artemis_ephemeris.json"
OUTPUT_FILE = "artemis_orion.dat"

LAUNCH_TIME = datetime(2026, 4, 1, 22, 35, 0, tzinfo=timezone.utc)
LAUNCH_TS = int(LAUNCH_TIME.timestamp())

with open(DATA_FILE) as f:
    rows = json.load(f)

with open(OUTPUT_FILE, "w") as f:
    # Header: number of records, launch unix timestamp
    f.write(f"{len(rows)} {LAUNCH_TS}\n")
    # Columns: timestamp  x  y  z  velocity_km_s  range_rate_km_s
    for r in rows:
        f.write(
            f"{r['timestamp']} "
            f"{r['x']:.12e} {r['y']:.12e} {r['z']:.12e} "
            f"{r['velocity_km_s']:.12e} {r['range_rate_km_s']:.12e}\n"
        )

print(f"Wrote {len(rows)} records to {OUTPUT_FILE}")
print(f"Launch timestamp: {LAUNCH_TS}")
