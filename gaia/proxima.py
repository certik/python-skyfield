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

if args.plot:
    import matplotlib.pyplot as plt

    # === 2. Plot 1: Proper-motion vector ===
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.arrow(0, 0, df['pmra'].iloc[0], df['pmdec'].iloc[0],
             head_width=0.1, head_length=0.2, fc='blue', ec='blue')
    ax.set_xlabel('Proper motion in RA (mas/yr)')
    ax.set_ylabel('Proper motion in Dec (mas/yr)')
    ax.set_title('Proxima Centauri Proper Motion Vector')
    ax.grid(True)
    ax.axis('equal')
    plt.show()

    # === 3. Plot 2: Projected sky path over next 100 years ===
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
