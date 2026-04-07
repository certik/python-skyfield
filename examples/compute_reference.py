#!/usr/bin/env python3
"""
Compute Sun/Moon reference alt/az/dist/diam using Skyfield + DE441s.

Uses the same inputs as our Fortran pipeline:
  - DE441s ephemeris
  - JPL EOP2 polar motion and delta_T
  - IAU 2006/2000A precession-nutation (Skyfield's default)

The resulting values serve as the reference for test_de441_horizons.f90
and test_compute_from_scratch.f90.  Both Skyfield (equinox-based) and
our Fortran (CIO-based) implement the same IAU standard, so they should
agree to sub-milliarcsecond level.

Test cases:
  1. 40°N, 0°E, 0 m — 2025-01-01 12:00 UTC
  2. Fredericksburg TX (30.275°N, 98.872°W, 556 m) — 2024-04-08 18:35:07 UTC
"""
import numpy as np
from skyfield.api import load, wgs84


SOLAR_RADIUS_KM = 695700.0
MOON_RADIUS_KM  = 1737.4
AU_KM           = 149597870.7


def load_eop2(filename, target_mjd_utc, leap_seconds):
    """Read polar motion and delta_T from JPL EOP2 for a given UTC date.

    Returns (xp_as, yp_as, delta_t_s, dX_mas, dY_mas).
    """
    mjd_tai = target_mjd_utc + leap_seconds / 86400.0
    in_data = False
    prev = None
    with open(filename) as f:
        for line in f:
            if line.strip().startswith('EOP2='):
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.strip().split(',')
            if len(parts) < 7:
                continue
            mjd = float(parts[0])
            if mjd > mjd_tai:
                if prev is None:
                    break
                frac = (mjd_tai - float(prev[0])) / (mjd - float(prev[0]))
                vals = []
                for i in (1, 2, 3, 5, 6):
                    vals.append(float(prev[i]) + frac *
                                (float(parts[i]) - float(prev[i])))
                xp_as = vals[0] / 1000.0
                yp_as = vals[1] / 1000.0
                delta_t = 32.184 + vals[2] / 1000.0
                dX_mas, dY_mas = vals[3], vals[4]
                return xp_as, yp_as, delta_t, dX_mas, dY_mas
            prev = parts
    raise RuntimeError("Date not found in EOP2 file")


def compute_reference(lat_deg, lon_deg, elev_m,
                      utc_year, utc_month, utc_day,
                      utc_hour, utc_minute, utc_second,
                      leap_seconds, eop_file, ephem_file):
    """Compute Sun/Moon alt/az/dist/diam using Skyfield with EOP2 inputs."""
    mjd_utc = (julian_day(utc_year, utc_month, utc_day) - 2400001
               + (utc_hour * 3600 + utc_minute * 60 + utc_second) / 86400.0)
    xp_as, yp_as, delta_t, _, _ = load_eop2(eop_file, mjd_utc, leap_seconds)

    jd_utc = julian_day(utc_year, utc_month, utc_day) + \
        (utc_hour * 3600 + utc_minute * 60 + utc_second) / 86400.0 - 0.5

    ts = load.timescale()
    ts.delta_t_function = lambda tt: delta_t
    ts.polar_motion_table = (
        np.array([jd_utc - 1.0, jd_utc + 1.0]),
        np.array([xp_as, xp_as]),
        np.array([yp_as, yp_as]),
    )

    t = ts.utc(utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_second)
    eph = load(ephem_file)
    earth, sun, moon = eph['earth'], eph['sun'], eph['moon']
    obs = earth + wgs84.latlon(lat_deg, lon_deg, elev_m)

    app_sun  = obs.at(t).observe(sun).apparent()
    app_moon = obs.at(t).observe(moon).apparent()

    alt_s, az_s, dist_s = app_sun.altaz()
    alt_m, az_m, dist_m = app_moon.altaz()

    sun_diam  = 2 * np.degrees(np.arcsin(
        SOLAR_RADIUS_KM / (dist_s.au * AU_KM))) * 3600
    moon_diam = 2 * np.degrees(np.arcsin(
        MOON_RADIUS_KM / (dist_m.au * AU_KM))) * 3600

    return {
        'delta_t': delta_t,
        'xp_as': xp_as, 'yp_as': yp_as,
        'sun_alt': alt_s.degrees, 'sun_az': az_s.degrees,
        'sun_dist': dist_s.au, 'sun_diam': sun_diam,
        'moon_alt': alt_m.degrees, 'moon_az': az_m.degrees,
        'moon_dist': dist_m.au, 'moon_diam': moon_diam,
    }


def julian_day(year, month, day):
    """Compute Julian Day Number for a calendar date."""
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + A // 4
    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524


def main():
    cases = [
        {
            'label': 'Test 1: 40°N, 0°E, 0 m — 2025-01-01 12:00:00 UTC',
            'lat_deg': 40.0, 'lon_deg': 0.0, 'elev_m': 0.0,
            'utc_year': 2025, 'utc_month': 1, 'utc_day': 1,
            'utc_hour': 12, 'utc_minute': 0, 'utc_second': 0,
            'leap_seconds': 37,
        },
        {
            'label': 'Test 2: Fredericksburg TX — 2024-04-08 18:35:07 UTC (eclipse)',
            'lat_deg': 30.2752011, 'lon_deg': -98.8719843, 'elev_m': 556.0,
            'utc_year': 2024, 'utc_month': 4, 'utc_day': 8,
            'utc_hour': 18, 'utc_minute': 35, 'utc_second': 7,
            'leap_seconds': 37,
        },
    ]

    for case in cases:
        label = case.pop('label')
        ref = compute_reference(**case, eop_file='latest_eop2.long',
                                ephem_file='de441s.bsp')

        print(f"Skyfield reference: {label}")
        print(f"delta_T = {ref['delta_t']:.6f} s")
        print(f"PM: xp = {ref['xp_as']:.6f}\", yp = {ref['yp_as']:.6f}\"")
        print()
        print(f"Sun  alt  = {ref['sun_alt']:.15f}°")
        print(f"Sun  az   = {ref['sun_az']:.15f}°")
        print(f"Sun  dist = {ref['sun_dist']:.17e} AU")
        print(f"Sun  diam = {ref['sun_diam']:.15f}\"")
        print(f"Moon alt  = {ref['moon_alt']:.15f}°")
        print(f"Moon az   = {ref['moon_az']:.15f}°")
        print(f"Moon dist = {ref['moon_dist']:.17e} AU")
        print(f"Moon diam = {ref['moon_diam']:.15f}\"")

        print()
        print("Fortran-ready constants:")
        for name, key, fmt in [
            ('ref_sun_alt',   'sun_alt',   '.15f'),
            ('ref_sun_az',    'sun_az',    '.15f'),
            ('ref_sun_dist',  'sun_dist',  '.17e'),
            ('ref_sun_diam',  'sun_diam',  '.15f'),
            ('ref_moon_alt',  'moon_alt',  '.15f'),
            ('ref_moon_az',   'moon_az',   '.15f'),
            ('ref_moon_dist', 'moon_dist', '.17e'),
            ('ref_moon_diam', 'moon_diam', '.15f'),
        ]:
            val = format(ref[key], fmt)
            print(f"real(dp), parameter :: {name:15s} = {val}_dp")
        print()
        print('=' * 70)
        print()


if __name__ == '__main__':
    main()
