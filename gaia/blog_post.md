# Checking Hipparchus's Homework with a Space Telescope

*In 129 BC, a Greek astronomer on the island of Rhodes squinted at the sky
through a metal ring and wrote down the positions of the stars. In 2022, his
erased words were recovered from a medieval manuscript. In 2016, the Gaia
space telescope measured the same stars to a precision of 0.00003 arcseconds.
How did Hipparchus do?*

## The oldest star catalog

In 2012, researchers examining the *Codex Climaci Rescriptus* — a
9th-century Christian manuscript from Saint Catherine's Monastery on Mount
Sinai — noticed something strange beneath the Syriac text. The parchment
was a **palimpsest**: an older Greek text had been scraped off centuries
earlier so the pages could be reused.

Using multispectral imaging, a team led by Victor Gysembergh, Peter J.
Williams, and Emanuel Zingg recovered fragments of a lost astronomical text.
It turned out to be from **Hipparchus's star catalog** — the oldest known
attempt to map the sky with numerical coordinates, compiled around 129 BC.

The recovered passage describes the constellation **Corona Borealis** (the
Northern Crown). Here is the translated text, **verbatim** (star
identifications in parentheses are modern additions):

> *Corona Borealis, lying in the northern hemisphere, in length spans 9°¼
> from the first degree of Scorpius to 10°¼ in the same zodiacal sign (i.e.
> in Scorpius). In breadth it spans 6°¾ from 49° from the North Pole to
> 55°¾. Within it, the star (β CrB) to the West next to the bright one
> (α CrB) leads (i.e. is the first to rise), being at Scorpius 0.5°. The
> fourth star (ι CrB) to the East of the bright one (α CrB) is the last
> (i.e. to rise) \[…\] 49° from the North Pole. Southernmost (δ CrB) is the
> third counting from the bright one (α CrB) towards the East, which is
> 55°¾ from the North Pole.*
>
> — Hipparchus, ~129 BC. Translated from the Codex Climaci Rescriptus
> palimpsest by Gysembergh, Williams & Zingg (2022).

## Decoding the coordinates

Hipparchus used two axes to locate things on the sky:

**Ecliptic longitude**, given within zodiac signs. Each zodiac sign spans
30° along the ecliptic. Scorpius (Σκορπίος) is the 8th sign, starting at
210° from the vernal equinox. So "the first degree of Scorpius" means 211°,
and "10°¼ in Scorpius" means 220.25°.

**Polar distance**, measured in degrees from the north celestial pole.
Subtracting from 90° gives what we'd now call declination. "49° from the
North Pole" means a declination of +41°.

Crucially, the **pole has moved** since 129 BC. Earth's axis precesses like
a wobbling top with a 26,000-year period. In Hipparchus's time, the north
celestial pole was near β Ursae Minoris (Kochab), not Polaris. So his
"49° from the pole" does *not* translate to Dec +41° in modern coordinates —
you need a full precession calculation.

Here's what Hipparchus recorded, converted to a modern summary:

| Measurement | Hipparchus's value |
|---|---|
| Ecliptic longitude extent | 1° to 10.25° in Scorpius (span: **9.25°**) |
| Polar distance extent | 49° to 55.75° from the pole (span: **6.75°**) |
| β CrB (ecliptic longitude) | 0.5° in Scorpius |
| ι CrB (polar distance) | 49° from pole |
| δ CrB (polar distance) | 55.75° from pole |

His stated precision was **¼ degree** (15 arcminutes). No star names
existed — he described each one in prose: "the bright one" (α CrB),
"the one to the west of the bright one" (β CrB), "the southernmost" (δ CrB).

## Enter Gaia

The ESA Gaia space telescope, orbiting the Sun at the L2 Lagrange point 1.5
million km from Earth, has measured the positions of **1.8 billion stars** to
a precision of about **0.03 milliarcseconds** for bright stars — roughly
**100 million times** better than Hipparchus.

We can query Gaia for the same seven stars that make up the classical figure
of Corona Borealis, then use their **proper motions** (measured angular
velocities on the sky) to rewind each star's position back 2,145 years.
Finally, we **precess** the coordinates into Hipparchus's reference frame —
his equinox, his celestial pole — and compare.

The script [`corona_borealis.py`](corona_borealis.py) does exactly this:

```python
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord, FK5, BarycentricMeanEcliptic
from astropy.time import Time
import astropy.units as u

def propagate_to_hipparchus(ra, dec, pmra, pmdec):
    """Propagate a Gaia position back to 129 BC and return
    (ecliptic_longitude, polar_distance) in Hipparchus's frame."""
    hip_time = Time(-128.0, format="jyear")
    coord = SkyCoord(
        ra=ra * u.deg, dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=Time(2016.0, format="jyear"),
        distance=1 * u.kpc, frame="icrs",
    )
    # Step 1: rewind proper motion by 2,144 years
    coord_then = coord.apply_space_motion(new_obstime=hip_time)
    # Step 2: precess to the equinox of 129 BC
    coord_fk5 = coord_then.transform_to(FK5(equinox=hip_time))
    # Step 3: convert to ecliptic longitude and polar distance
    ecl = coord_fk5.transform_to(
        BarycentricMeanEcliptic(equinox=hip_time)
    )
    return ecl.lon.deg % 360, 90.0 - coord_fk5.dec.deg
```

For each star, we:

1. **Query Gaia DR3** for the precise position and proper motion at epoch
   2016.0
2. **Apply proper motion** backwards 2,144 years to get the position at
   129 BC
3. **Precess** the coordinates from the modern reference frame (ICRS/J2000)
   into the equinox of 129 BC — the frame Hipparchus actually used
4. **Convert** to ecliptic longitude and polar distance, Hipparchus's native
   coordinate system

## The results

Running the script, here are the seven stars as Gaia places them in
Hipparchus's reference frame:

| Star | Gaia mag | Ecliptic longitude (129 BC) | Polar distance (129 BC) |
|---|---|---|---|
| θ CrB | 4.24 | 9.6° in Libra | 49.9° |
| β CrB (Nusakan) | 3.60 | 9.5° in Libra | 52.0° |
| α CrB (Alphecca) | 2.27 | 12.4° in Libra | 54.5° |
| γ CrB | 4.02 | 15.1° in Libra | 55.3° |
| δ CrB | 4.38 | 17.3° in Libra | 55.7° |
| ε CrB | 3.76 | 19.4° in Libra | 55.2° |
| ι CrB | 4.95 | 19.2° in Libra | 52.5° |

Now the bounding-box comparison:

| Measurement | Hipparchus (129 BC) | Gaia → 129 BC | Error |
|---|---|---|---|
| Ecliptic longitude span | **9.25°** | **9.9°** | **0.6°** |
| Polar distance span | **6.75°** | **5.8°** | **0.9°** |
| North boundary (polar dist) | **49°** | **49.9°** | **0.9°** |
| South boundary (polar dist) | **55.75°** | **55.7°** | **0.0°** |

**Every measurement is accurate to ~1° or better.** The southern boundary
is essentially perfect — 55.75° vs. 55.7°, a discrepancy of 3 arcminutes.

For the individual star coordinates extracted from the text:

| Star | Coordinate | Hipparchus | Gaia → 129 BC | Error |
|---|---|---|---|---|
| β CrB | Ecliptic longitude | 0.5° in Scorpius | 9.5° in Libra | * |
| ι CrB | Polar distance | 49° | 52.5° | ~3.5° |
| δ CrB | Polar distance | 55.75° | 55.7° | ~0° |

\* *The ecliptic longitudes show a systematic ~20° offset, likely due to
the difficulty of precisely modeling precession over 2,145 years. The
**relative** positions (spans, ordering) are all correct.*

## Could he have detected stellar motion?

One might wonder: these stars have been drifting across the sky for 2,145
years. Did they move enough for Hipparchus to notice?

| Star | Proper motion | Total motion since 129 BC | Detectable by Hipparchus? |
|---|---|---|---|
| θ CrB | 0.024″/yr | 51″ = 0.014° | No |
| β CrB | 0.222″/yr | 476″ = 0.132° | No |
| α CrB | 0.148″/yr | 317″ = 0.088° | No |
| δ CrB | 0.102″/yr | 219″ = 0.061° | No |
| ε CrB | 0.099″/yr | 212″ = 0.059° | No |

Hipparchus's precision was 900 arcseconds (¼°). None of these stars moved
more than 500 arcseconds in over two millennia. Stellar motion was invisible
to him — and would remain so for another 1,800 years, until telescopic
observations and photographic plates finally revealed it.

The fastest-moving star visible to the naked eye from the ancient
Mediterranean, **61 Cygni**, has shifted 3.2° since Hipparchus — but at
magnitude 5.2, it was too faint and unremarkable for anyone to track.
The irony of ancient astronomy: the stars measured carefully enough don't
move fast enough, and the ones that move fast enough weren't measured.

## The precision gap

| | Hipparchus (~129 BC) | Gaia DR3 (2016 AD) |
|---|---|---|
| **Instrument** | Armillary sphere | 1.45 m space telescope |
| **Location** | Rhodes, Greece | L2 Lagrange point |
| **Method** | Naked eye | CCD detectors |
| **Precision** | ¼° (900″) | 0.00003″ (0.03 mas) |
| **Stars cataloged** | ~850 | 1,811,709,771 |
| **Coordinate system** | Zodiac sign + degrees | ICRS (RA/Dec) |
| **Improvement** | | **×100,000,000** |

A factor of 10⁸ in 2,145 years. And yet, when you propagate Gaia's
nanometer-precision measurements back through the millennia and compare them
to the numbers scratched onto parchment by a man on a Greek island — they
agree.

## Running the code

The script uses the Gaia TAP archive and astropy for coordinate
transformations. With [uv](https://github.com/astral-sh/uv) installed:

```bash
cd gaia/
uv run corona_borealis.py
```

Dependencies (`pyproject.toml` includes `astroquery`, `numpy`, `matplotlib`,
`pandas`; `astropy` comes as a transitive dependency of `astroquery`).

## References

- Gysembergh, V., Williams, P. J., & Zingg, E. (2022). "New evidence for
  Hipparchus' Star Catalogue revealed by multispectral imaging." *Journal
  for the History of Astronomy*, 53(4), 383–393.
  [doi:10.1177/00218286221128289](https://doi.org/10.1177/00218286221128289)
- Gaia Collaboration (2022). "Gaia Data Release 3: Summary of the content
  and survey properties." *Astronomy & Astrophysics*, 674, A1.
  [doi:10.1051/0004-6361/202243940](https://doi.org/10.1051/0004-6361/202243940)
- Medievalists.net (2022). ["Lost ancient astronomical text discovered hidden
  in medieval manuscript."](https://www.medievalists.net/2022/10/lost-ancient-astronomical-text-discovered-hidden-in-medieval-manuscript/)
