# What Can You Measure from a Backyard?

## Setup

A single observer at a known latitude and longitude measures the **altitude and
azimuth** of the Sun and Moon over one year, with **1 arcminute (60″) Gaussian
noise** per measurement, roughly every clear day/night — about 3000 Sun and 3000
Moon observations total.

We assume two length scales are known from independent measurements:

- **R_E** — the Earth's radius (6 378 km)
- **a_E** — the Earth–Sun distance, i.e. the astronomical unit (1.496 × 10⁸ km)

We also assume the Moon's gravitational parameter **GM_Moon** is known (see
below for why).

Everything else — orbital elements, masses of the Sun and Earth, the Moon's
distance — we attempt to derive from the observations alone, using an Extended
Kalman Filter with a 3-body (Sun–Earth–Moon) N-body integrator.

---

## What the observations constrain well

### 1. All orbital angles (both orbits)

The angular positions of the Sun and Moon on the sky directly determine:

| Parameter | Body  | Typical precision |
|-----------|-------|-------------------|
| e (eccentricity) | Earth | < 0.1% |
| i (inclination)  | Earth | < 0.01° |
| Ω (longitude of ascending node) | Earth | < 0.1° |
| ω (argument of perihelion) | Earth | < 0.1° |
| M₀ (mean anomaly at epoch) | Earth | < 0.1° |
| e, i, Ω, ω, M₀ | Moon | comparable |

These are the "shape and orientation" of each orbit. Angular measurements
determine angles very well.

### 2. Orbital periods (mean motions)

From the angular rate of the Sun and Moon across the sky:

- **T_E ≈ 365.25 days** (Earth orbital period) — determined to ~0.01%
- **T_M ≈ 27.32 days** (Moon orbital period) — determined to ~0.01%

### 3. Moon distance (a_M) — via diurnal parallax

As the Earth rotates, the observer moves by up to R_E × cos(lat) ≈ 4 900 km.
This produces a **diurnal parallax** of the Moon:

    Moon parallax ≈ R_E / a_M ≈ 57 arcminutes

This is ~57× larger than the measurement noise. With thousands of observations
the Moon's distance is determined to high precision:

    a_M ≈ 384 400 km (precision < 0.1%)

### 4. μ_EM = GM_Earth + GM_Moon — from Kepler's third law

Knowing both a_M and T_M, Kepler's third law gives:

    μ_EM = (2π)² × a_M³ / T_M² ≈ 4.086 × 10⁵ km³/s²

This is determined as well as a_M (< 0.1%).

### 5. GM_Earth — by subtracting GM_Moon

Since GM_Moon is assumed known (4 903 km³/s²):

    GM_Earth = μ_EM − GM_Moon ≈ 3.986 × 10⁵ km³/s²

### 6. μ_SE = GM_Sun + GM_Earth — from Kepler's third law (Earth orbit)

Knowing a_E (assumed) and T_E (measured):

    μ_SE = (2π)² × a_E³ / T_E²

Since T_E is very precisely measured and a_E is assumed known, μ_SE is
well-determined. Then:

    GM_Sun = μ_SE − GM_Earth ≈ 1.327 × 10¹¹ km³/s²

---

## What the observations cannot constrain

### 1. The astronomical unit (a_E)

This is the central difficulty we encountered. **Pure angular observations from a
single location cannot determine the absolute Sun–Earth distance.**

From angular measurements, the Earth's mean motion n_E = 2π/T_E is known
precisely. Kepler's third law then gives:

    n_E² = μ_SE / a_E³

This is one equation in two unknowns (μ_SE and a_E). Any pair satisfying this
relation produces **identical** angular motion of the Sun. The system is scale-
invariant: multiplying all distances and masses by constants that preserve n_E
leaves all angles unchanged.

**Can the Moon break this degeneracy?** Almost, but not quite:

- The dominant solar perturbations on the Moon (evection, variation, annual
  equation) arise from the **quadrupole** (P₂) term of the Sun's tidal
  potential. They scale as μ_SE/a_E³ — the same combination as n_E². So they
  add no new information about a_E.

- The **parallactic inequality** (~125″ amplitude, synodic period) arises from
  the **octupole** (P₃) term and scales as μ_SE/a_E⁴ — one extra power of
  1/a_E. In principle, comparing its amplitude to the variation gives the ratio
  a_M/a_E, hence a_E.

- In practice, extracting a_E this way requires measuring a ~0.3% difference in
  how μ_SE and a_E affect the Jacobian. At 60″ noise this yields ~1–2%
  precision on a_E at best — and the EKF's finite-difference Jacobian struggles
  to capture this octupole-level sensitivity numerically.

**Historical context:** Flamsteed (1670s) and Tobias Mayer (1750s) attempted
this approach. Even Mayer, with ~15″ lunar observations, could only constrain
the solar parallax to ~10″ ± 3″ (true: 8.8″). The AU was not precisely measured
until radar ranging to Venus in the 1960s.

### 2. Individual Moon mass (GM_Moon)

The Moon's orbit depends on μ_EM = GM_Earth + GM_Moon, not on each mass
individually. The one signal that separates them — the monthly wobble of the
Sun's apparent position due to the Earth–Moon barycenter offset — has amplitude:

    δθ ≈ (a_M × GM_Moon / μ_EM) / a_E ≈ 6.2″

This is 10× below the single-measurement noise. Statistically it might reach
~5σ with 3000 Sun observations, but is heavily degenerate with Earth orbital
elements. At 60″ noise, GM_Moon must be assumed known.

### 3. Solar parallax directly

The solar parallax (topocentric shift of the Sun due to R_E):

    Solar parallax = R_E / a_E ≈ 8.8″

This is ~7× below the measurement noise and undetectable.

---

## Summary table

| Quantity | Determined? | How | Precision |
|----------|------------|-----|-----------|
| Earth orbital angles (e, i, Ω, ω, M₀) | ✅ Yes | Sun alt/az over 1 year | < 0.1° |
| Moon orbital angles (e, i, Ω, ω, M₀) | ✅ Yes | Moon alt/az over 1 year | < 0.1° |
| Earth orbital period (T_E) | ✅ Yes | Sun angular rate | < 0.01% |
| Moon orbital period (T_M) | ✅ Yes | Moon angular rate | < 0.01% |
| Moon distance (a_M) | ✅ Yes | Diurnal parallax (57′ ≫ 60″) | < 0.1% |
| μ_EM = GM_Earth + GM_Moon | ✅ Yes | Kepler: a_M + T_M | < 0.1% |
| GM_Earth | ✅ Yes | μ_EM − GM_Moon (assumed) | < 0.1% |
| GM_Sun | ✅ Yes | μ_SE − GM_Earth (requires a_E assumed) | < 0.1% |
| Earth–Sun distance (a_E = AU) | ❌ No | Would need solar parallax (8.8″) or precise PI extraction | ~1–2% at best |
| GM_Moon individually | ❌ No | EMB wobble signal is only 6.2″ | unresolvable |
| Solar parallax | ❌ No | 8.8″ signal in 60″ noise | unresolvable |

---

## What must be assumed

1. **R_E (Earth radius)** — the only "ruler" connecting angular measurements to
   absolute distances. Provides the Moon distance through parallax.

2. **a_E (the AU)** — cannot be derived from 1-arcminute angular observations.
   Must be taken from an independent measurement (historically: transit of
   Venus, asteroid parallax, or radar ranging).

3. **GM_Moon** — the Earth/Moon mass split cannot be resolved at this noise
   level. The total μ_EM is well-determined, but splitting it requires ~1″
   measurements.

## What can be derived

Given those three assumed inputs, backyard 1-arcminute observations over one
year are sufficient to determine:

- The complete orbital geometry of both the Earth and Moon
- The orbital periods
- The distance to the Moon
- The mass of the Earth
- The mass of the Sun
- The lunar perturbation structure (evection, variation, annual equation,
  parallactic inequality) — all detected as the N-body model naturally
  reproduces them

This is a remarkable result: **three assumed constants and a year of amateur
observations yield the masses and orbits of the entire Sun–Earth–Moon system.**
