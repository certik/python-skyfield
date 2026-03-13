import argparse

from astroquery.gaia import Gaia
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Query Gaia DR3 for Proxima Centauri")
parser.add_argument("--plot", action="store_true", help="Show plots")
args = parser.parse_args()

# === 1. Download all data from gaiadr3.gaia_source ===
query = """
SELECT *
FROM gaiadr3.gaia_source
WHERE source_id = 5853498713190525696
"""

print("Querying Gaia Archive for Proxima Centauri...")
job = Gaia.launch_job(query)
results = job.get_results()

df = results.to_pandas()

# Print every column as a name/value table
print("\n=== Gaia DR3 Data for Proxima Centauri ===\n")
row = df.iloc[0]
max_name = max(len(str(c)) for c in row.index)
for col, val in row.items():
    print(f"  {col:<{max_name}}  {val}")

# Quick distance in light-years
parallax_mas = df['parallax'].iloc[0]
distance_pc = 1000 / parallax_mas
distance_ly = distance_pc * 3.26156
print(f"\nDistance: {distance_pc:.4f} pc ≈ {distance_ly:.4f} light-years")
print(f"Parallax uncertainty: ±{df['parallax_error'].iloc[0]:.4f} mas")

# === 2. RA/Dec as a function of time ===
ra0 = df['ra'].iloc[0]           # deg at epoch 2016.0
dec0 = df['dec'].iloc[0]         # deg at epoch 2016.0
pmra = df['pmra'].iloc[0]        # mas/yr (includes cos(dec) factor)
pmdec = df['pmdec'].iloc[0]      # mas/yr
parallax = df['parallax'].iloc[0]  # mas
rv = df['radial_velocity'].iloc[0]  # km/s

cos_dec = np.cos(np.radians(dec0))

# Propagate position yearly from 2014 to 2030, including parallax wobble.
# Earth's orbital phase (simplified circular orbit, vernal equinox ~ March 20)
ref_epoch = 2016.0
years = np.arange(2014, 2031, 0.25)  # quarterly steps

print("\n=== RA/Dec as a function of time ===\n")
print(f"  Reference position (epoch {ref_epoch}):")
print(f"    RA  = {ra0:.10f}°")
print(f"    Dec = {dec0:.10f}°\n")

# Ecliptic longitude of the star (approximate, for parallax wobble direction)
ecl_lon_rad = np.radians(df['ecl_lon'].iloc[0])
ecl_lat_rad = np.radians(df['ecl_lat'].iloc[0])

print(f"  {'Year':>8}  {'RA (°)':>18}  {'Dec (°)':>18}  {'ΔRA (arcsec)':>14}  {'ΔDec (arcsec)':>14}")
print(f"  {'':->8}  {'':->18}  {'':->18}  {'':->14}  {'':->14}")

for t in years:
    dt = t - ref_epoch

    # Proper motion contribution (mas → deg)
    dra_pm = pmra * dt / (3600_000 * cos_dec)
    ddec_pm = pmdec * dt / 3600_000

    # Parallax wobble (simplified: Earth at angular position 2π*(t mod 1))
    # The wobble in ecliptic coords is: Δlon = parallax/cos(lat) * sin(phase)
    # and Δlat = parallax * sin(lat) * cos(phase), but we compute in RA/Dec
    # using the standard approximation for the parallactic displacement.
    phase = 2 * np.pi * (t - 0.22)  # ~0.22 for vernal equinox offset (~March 20)
    # Parallactic factors (approximate)
    p_ra = parallax * (np.sin(phase) * np.sin(ecl_lon_rad)
           - np.cos(phase) * np.cos(ecl_lon_rad) * np.sin(ecl_lat_rad))
    p_dec = parallax * (np.sin(phase) * np.cos(ecl_lon_rad) * np.sin(np.radians(dec0))
            - np.cos(phase) * (np.sin(ecl_lat_rad) * np.cos(np.radians(dec0))
            - np.cos(ecl_lat_rad) * np.sin(np.radians(dec0)) * np.sin(ecl_lon_rad)))
    # Simple approximation: project parallax ellipse
    dra_plx = p_ra / (3600_000 * cos_dec)
    ddec_plx = p_dec / 3600_000

    ra_t = ra0 + dra_pm + dra_plx
    dec_t = dec0 + ddec_pm + ddec_plx

    dra_arcsec = (ra_t - ra0) * 3600 * cos_dec
    ddec_arcsec = (dec_t - dec0) * 3600

    print(f"  {t:8.2f}  {ra_t:18.10f}  {dec_t:18.10f}  {dra_arcsec:+14.4f}  {ddec_arcsec:+14.4f}")

# === 3. Radial distance ===
print("\n=== Radial distance ===\n")
distance_pc = 1000 / parallax
distance_ly = distance_pc * 3.26156
distance_au = distance_pc * 206265
print(f"  Parallax:  {parallax:.4f} ± {df['parallax_error'].iloc[0]:.4f} mas")
print(f"  Distance:  {distance_pc:.6f} pc")
print(f"           = {distance_ly:.6f} light-years")
print(f"           = {distance_au:.0f} AU")

# === 4. Radial velocity ===
print("\n=== Radial velocity ===\n")
print(f"  Radial velocity:  {rv:.4f} ± {df['radial_velocity_error'].iloc[0]:.4f} km/s")
print(f"  (negative = approaching us)")
print(f"  Method: Doppler shift of spectral lines (RVS spectrograph)")
print(f"  Transits used: {df['rv_nb_transits'].iloc[0]}")

# Tangential velocity for context
vt = 4.74047 * (df['pm'].iloc[0] / 1000) * distance_pc
vtot = np.sqrt(vt**2 + rv**2)
print(f"\n  Tangential velocity (from proper motion × distance): {vt:.2f} km/s")
print(f"  Total 3D space velocity: {vtot:.2f} km/s")

if args.plot:
    import matplotlib.pyplot as plt

    # === 3. Plot 1: Proper-motion vector ===
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.arrow(0, 0, df['pmra'].iloc[0], df['pmdec'].iloc[0],
             head_width=0.1, head_length=0.2, fc='blue', ec='blue')
    ax.set_xlabel('Proper motion in RA (mas/yr)')
    ax.set_ylabel('Proper motion in Dec (mas/yr)')
    ax.set_title('Proxima Centauri Proper Motion Vector')
    ax.grid(True)
    ax.axis('equal')
    plt.show()

    # === 4. Plot 2: Projected sky path over next 100 years ===
    years = np.linspace(0, 100, 101)
    pmra = df['pmra'].iloc[0]      # mas/yr
    pmdec = df['pmdec'].iloc[0]    # mas/yr

    ra_path = years * pmra / 1000   # arcsec
    dec_path = years * pmdec / 1000 # arcsec

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ra_path, dec_path, 'b-', label='Projected path')
    ax.scatter([0], [0], color='red', s=50, label='Position in 2016.0')
    ax.set_xlabel('Displacement in RA (arcsec)')
    ax.set_ylabel('Displacement in Dec (arcsec)')
    ax.set_title('Proxima Centauri Sky Path: Next 100 Years (proper motion only)')
    ax.grid(True)
    ax.legend()
    plt.show()
