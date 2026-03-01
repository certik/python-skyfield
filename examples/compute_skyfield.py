#!/usr/bin/env python3
"""
Compute Sun and Moon positions using Skyfield's public API only.

No hand-computation — just Skyfield calls. Computes altitude, azimuth,
and angular radius for the April 8, 2024 total solar eclipse as seen
from Fredericksburg, TX.
"""

from datetime import datetime
from math import asin, atan
from pytz import timezone
from skyfield.api import load, wgs84

SOLAR_RADIUS_KM = 696340.0
MOON_RADIUS_KM = 1737.1

eph = load('de440s.bsp')
earth, sun, moon = eph['earth'], eph['sun'], eph['moon']
ts = load.timescale()

observer = wgs84.latlon(30.274167, -98.871944, 516)
zone = timezone('US/Central')
time = zone.localize(datetime(2024, 4, 8, 13, 35, 10))
t = ts.from_datetime(time)

print("Location: Fredericksburg, TX")
print(f"Time:     {t.astimezone(zone)}")
print()

obs = (earth + observer).at(t)

# ── Sun ──
sun_app = obs.observe(sun).apparent()
sun_alt, sun_az, sun_dist = sun_app.altaz()
sun_radius_asin = asin(SOLAR_RADIUS_KM / sun_dist.km)
sun_radius_atan = atan(SOLAR_RADIUS_KM / sun_dist.km)

print("Sun:")
print(f"  Altitude: {sun_alt}")
print(f"  Azimuth:  {sun_az}")
print(f"  Distance: {sun_dist.km:.3f} km  ({sun_dist.au:.15f} AU)")
print(f"  Radius (asin): {sun_radius_asin:.16f} rad  ({sun_radius_asin * 180 / 3.141592653589793:.6f}°)")
print(f"  Radius (atan): {sun_radius_atan:.16f} rad  ({sun_radius_atan * 180 / 3.141592653589793:.6f}°)")
print()

# ── Moon ──
moon_app = obs.observe(moon).apparent()
moon_alt, moon_az, moon_dist = moon_app.altaz()
moon_radius_asin = asin(MOON_RADIUS_KM / moon_dist.km)
moon_radius_atan = atan(MOON_RADIUS_KM / moon_dist.km)

print("Moon:")
print(f"  Altitude: {moon_alt}")
print(f"  Azimuth:  {moon_az}")
print(f"  Distance: {moon_dist.km:.3f} km  ({moon_dist.au:.15f} AU)")
print(f"  Radius (asin): {moon_radius_asin:.16f} rad  ({moon_radius_asin * 180 / 3.141592653589793:.6f}°)")
print(f"  Radius (atan): {moon_radius_atan:.16f} rad  ({moon_radius_atan * 180 / 3.141592653589793:.6f}°)")
print()

# ── Compare with plot3.py reference values ──
ref = (67.31640064112162, 178.74688214701396, 0.0046479241546984506,
       67.31763737226208, 178.74983743582982, 0.004908042955826713)

print("=" * 60)
print("Comparison with plot3.py reference values")
print("=" * 60)
pairs = [
    ("Sun alt",    sun_alt.degrees,  ref[0]),
    ("Sun az",     sun_az.degrees,   ref[1]),
    ("Sun radius", sun_radius_asin,  ref[2]),
    ("Moon alt",   moon_alt.degrees, ref[3]),
    ("Moon az",    moon_az.degrees,  ref[4]),
    ("Moon radius",moon_radius_asin, ref[5]),
]
for label, computed, expected in pairs:
    diff = abs(computed - expected)
    ok = "✓" if diff < 1e-10 else "✗"
    print(f"  {ok} {label:12s}: computed={computed:.16f}  ref={expected:.16f}  diff={diff:.2e}")
