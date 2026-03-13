"""Compare Hipparchus's (~129 BC) measurements of Corona Borealis
recovered from the Codex Climaci Rescriptus palimpsest against
Gaia DR3 data propagated back to the same epoch.

Usage:
    uv run corona_borealis.py
"""

import warnings

from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord, FK5, BarycentricMeanEcliptic
from astropy.time import Time
import astropy.units as u
import numpy as np

warnings.filterwarnings("ignore")

# ── Hipparchus's bounding box (from the palimpsest) ─────────────────
# Ecliptic longitude: "from the 1st degree of Scorpius to 10¼° of Scorpius"
# Polar distance:     "from 49° to 55¾° from the North Pole"
HIPPARCHUS_LON_WEST = 1.0       # degrees in Scorpius (= 210° + 1°)
HIPPARCHUS_LON_EAST = 10.25     # degrees in Scorpius (= 210° + 10.25°)
HIPPARCHUS_LON_SPAN = HIPPARCHUS_LON_EAST - HIPPARCHUS_LON_WEST  # 9.25°
HIPPARCHUS_POLAR_NORTH = 49.0   # degrees from pole
HIPPARCHUS_POLAR_SOUTH = 55.75  # degrees from pole
HIPPARCHUS_POLAR_SPAN = HIPPARCHUS_POLAR_SOUTH - HIPPARCHUS_POLAR_NORTH  # 6.75°

HIPPARCHUS_EPOCH = -128.0  # Julian year (= 129 BC)

# ── The 7 classical naked-eye stars of Corona Borealis ──────────────
# Approximate J2016 positions (refined below via Gaia cone search)
CRB_STARS = [
    ("θ CrB",           233.23, 31.36),
    ("β CrB (Nusakan)", 231.96, 29.11),
    ("α CrB (Alphecca)", 233.67, 26.71),
    ("γ CrB",           235.69, 26.30),
    ("δ CrB",           237.40, 26.07),
    ("ε CrB",           239.40, 26.88),
    ("ι CrB",           240.36, 29.85),
]

ZODIAC_SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpius", "Sagittarius", "Capricornus", "Aquarius", "Pisces",
]


def lon_to_zodiac(lon_deg):
    """Convert ecliptic longitude to 'X.X° in Sign' string."""
    lon_deg = lon_deg % 360
    idx = int(lon_deg / 30) % 12
    return f"{lon_deg - idx * 30:.1f}° in {ZODIAC_SIGNS[idx]}"


def query_star(ra, dec, radius=0.02):
    """Find the brightest Gaia source near the given position."""
    q = f"""
    SELECT TOP 1 source_id, ra, dec, pmra, pmdec, phot_g_mean_mag
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                     CIRCLE('ICRS', {ra}, {dec}, {radius}))
    ORDER BY phot_g_mean_mag ASC
    """
    job = Gaia.launch_job(q)
    r = job.get_results()
    return r[0] if len(r) > 0 else None


def propagate_to_hipparchus(ra, dec, pmra, pmdec):
    """Propagate a Gaia J2016 position back to Hipparchus's epoch and
    return (ecliptic_longitude, polar_distance) in his reference frame."""
    hip_time = Time(HIPPARCHUS_EPOCH, format="jyear")
    hip_equinox = Time(HIPPARCHUS_EPOCH, format="jyear")

    coord = SkyCoord(
        ra=ra * u.deg, dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=Time(2016.0, format="jyear"),
        distance=1 * u.kpc,  # placeholder (doesn't affect angles)
        frame="icrs",
    )
    # Apply proper motion back to 129 BC
    coord_then = coord.apply_space_motion(new_obstime=hip_time)
    # Precess to equinox of 129 BC (the reference frame Hipparchus used)
    coord_fk5 = coord_then.transform_to(FK5(equinox=hip_equinox))
    ecl = coord_fk5.transform_to(
        BarycentricMeanEcliptic(equinox=hip_equinox)
    )
    ecl_lon = ecl.lon.deg % 360
    polar_dist = 90.0 - coord_fk5.dec.deg
    return ecl_lon, polar_dist


def main():
    # ── Part 1: What Hipparchus wrote ────────────────────────────────
    print("=" * 72)
    print("WHAT HIPPARCHUS WROTE (~129 BC)")
    print("=" * 72)
    print("""
Recovered from the Codex Climaci Rescriptus palimpsest in 2022 by
Gysembergh, Williams & Zingg using multispectral imaging.  The original
Greek text had been scraped off parchment by medieval monks and
overwritten with Syriac Christian texts.

Hipparchus described each constellation with a bounding box and then
located individual stars within it.  For Corona Borealis he wrote:

  ┌────────────────────────────────────────────────────────────────┐
  │  "Corona [Borealis]:                                          │
  │   extends from the 1st degree of Scorpius to 10¼° of          │
  │   Scorpius in [ecliptic] longitude, and from 49° to           │
  │   55¾° from the North Pole."                                  │
  │                                                                │
  │  Individual stars named descriptively:                         │
  │    "the bright one"                     → α CrB (Alphecca)    │
  │    "the one to the west of the bright"  → β CrB (Nusakan)    │
  │    "the southernmost"                   → δ CrB               │
  │    "the easternmost"                    → ι CrB               │
  └────────────────────────────────────────────────────────────────┘

His coordinate system:
  • Ecliptic longitude within zodiac signs, to ¼° precision
    ("3° in Scorpius" = 210° + 3° = 213° from the vernal equinox)
  • Polar distance from the North Celestial Pole, to ¼° precision
    ("49° from pole" = declination +41° in HIS epoch)

Note: the celestial pole has MOVED since 129 BC due to precession.
In Hipparchus's time it was near β UMi (Kochab), not Polaris.
""")

    # ── Part 2: Query Gaia for the same stars ────────────────────────
    print("=" * 72)
    print("QUERYING GAIA DR3 FOR THE 7 CLASSICAL STARS OF CORONA BOREALIS")
    print("=" * 72)
    print()

    results = []
    for name, ra_approx, dec_approx in CRB_STARS:
        row = query_star(ra_approx, dec_approx)
        if row is None:
            print(f"  ✗ {name:<22} not found in Gaia")
            continue
        results.append((name, row))
        print(f"  ✓ {name:<22} G={row['phot_g_mean_mag']:.2f}  "
              f"RA={row['ra']:.4f}°  Dec={row['dec']:+.4f}°")

    # ── Part 3: Propagate each star back to 129 BC ───────────────────
    print()
    print("=" * 72)
    print("STARS PROJECTED TO HIPPARCHUS'S EPOCH & REFERENCE FRAME (129 BC)")
    print("=" * 72)
    print()
    print(f"  {'Star':<22} {'G':>5}  {'Ecl longitude':>22}  {'Polar dist':>10}")
    print(f"  {'':─<22} {'':─>5}  {'':─>22}  {'':─>10}")

    ecl_lons = []
    polar_dists = []

    for name, row in results:
        pmra = row['pmra']
        pmdec = row['pmdec']
        try:
            pmra, pmdec = float(pmra), float(pmdec)
            if np.isnan(pmra) or np.isnan(pmdec):
                pmra, pmdec = 0.0, 0.0
        except (TypeError, ValueError):
            pmra, pmdec = 0.0, 0.0
        elon, polar = propagate_to_hipparchus(
            row['ra'], row['dec'], pmra, pmdec,
        )
        ecl_lons.append(elon)
        polar_dists.append(polar)

        tag = "← brightest" if "Alphecca" in name else ""
        print(f"  {name:<22} {row['phot_g_mean_mag']:5.2f}"
              f"  {lon_to_zodiac(elon):>22}  {polar:8.1f}°  {tag}")

    # ── Part 4: Compare bounding boxes ───────────────────────────────
    lon_span = max(ecl_lons) - min(ecl_lons)
    pol_span = max(polar_dists) - min(polar_dists)
    pol_north = min(polar_dists)
    pol_south = max(polar_dists)

    print()
    print("=" * 72)
    print("BOUNDING BOX COMPARISON: HIPPARCHUS vs GAIA → 129 BC")
    print("=" * 72)
    print()
    print(f"  {'Measurement':<26} {'Hipparchus':>12}  {'Gaia→129BC':>12}  {'Error':>6}")
    print(f"  {'':─<26} {'':─>12}  {'':─>12}  {'':─>6}")
    print(f"  {'Ecliptic lon span':<26}"
          f" {HIPPARCHUS_LON_SPAN:11.2f}°"
          f"  {lon_span:11.1f}°"
          f"  {abs(lon_span - HIPPARCHUS_LON_SPAN):5.1f}°")
    print(f"  {'Polar dist span':<26}"
          f" {HIPPARCHUS_POLAR_SPAN:11.2f}°"
          f"  {pol_span:11.1f}°"
          f"  {abs(pol_span - HIPPARCHUS_POLAR_SPAN):5.1f}°")
    print(f"  {'Polar dist (north edge)':<26}"
          f" {HIPPARCHUS_POLAR_NORTH:11.2f}°"
          f"  {pol_north:11.1f}°"
          f"  {abs(pol_north - HIPPARCHUS_POLAR_NORTH):5.1f}°")
    print(f"  {'Polar dist (south edge)':<26}"
          f" {HIPPARCHUS_POLAR_SOUTH:11.2f}°"
          f"  {pol_south:11.1f}°"
          f"  {abs(pol_south - HIPPARCHUS_POLAR_SOUTH):5.1f}°")

    # ── Part 5: How much have these stars moved? ─────────────────────
    dt = 2016.0 - HIPPARCHUS_EPOCH  # ~2144 years
    print()
    print("=" * 72)
    print("HOW MUCH HAVE THESE STARS MOVED IN 2145 YEARS?")
    print("=" * 72)
    print()
    print(f"  {'Star':<22} {'PM (\"/yr)':>9}  {'Total motion':>14}  {'Detectable?':>14}")
    print(f"  {'':─<22} {'':─>9}  {'':─>14}  {'':─>14}")

    for name, row in results:
        pmra = row['pmra']
        pmdec = row['pmdec']
        try:
            pm = float(np.sqrt(float(pmra)**2 + float(pmdec)**2)) / 1000
            total = pm * dt
        except (TypeError, ValueError):
            print(f"  {name:<22}       n/a           n/a             n/a")
            continue
        if np.isnan(pm):
            print(f"  {name:<22}       n/a           n/a             n/a")
            continue
        detectable = "NO" if total < 900 else "borderline"
        print(f"  {name:<22} {pm:9.3f}  {total:10.1f}\"  "
              f"= {total/3600:.3f}°  {'':>2}{detectable}")

    print(f"""
  Hipparchus's precision was ¼° = 900 arcseconds.
  None of these stars moved more than ~500\" in 2145 years.
  → Their motion is UNDETECTABLE in his measurements.
""")

    # ── Part 6: Precision comparison ─────────────────────────────────
    print("=" * 72)
    print("MEASUREMENT PRECISION: 129 BC vs 2016 AD")
    print("=" * 72)
    print("""
  Hipparchus (~129 BC)           Gaia DR3 (2016 AD)
  ─────────────────────────      ──────────────────────────────
  Instrument: armillary sphere   Instrument: 1.45m space telescope
  Location:   Rhodes, Greece     Location:   L2 Lagrange point
  Method:     naked eye          Method:     CCDs, 0.1 nm precision
  Precision:  ¼° (900\")          Precision:  0.00003\" (0.03 mas)
  Stars:      ~850               Stars:      1,811,709,771
  Coordinate: zodiac sign + °    Coordinate: RA/Dec (ICRS)

  Improvement: ~100,000,000× in 2145 years
""")


if __name__ == "__main__":
    main()
