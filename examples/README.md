# Skyfield Plot Examples

Astronomy plots using [Skyfield](https://rhodesmill.org/skyfield/) and
[Matplotlib](https://matplotlib.org/).

## Prerequisites

- [uv](https://docs.astral.sh/uv/)

That's it — `uv` handles everything else (Python, dependencies, virtual env).

## Quick Start

```bash
cd examples/
uv run plot.py
```

On first run, `uv` will create a virtual environment, install all dependencies
(including the local Skyfield from the parent directory), and execute the script.
Skyfield will also auto-download the `de421.bsp` ephemeris file (~17 MB) on
first use.

## Scripts

### `plot.py` — Solar System Orbits

Plots the orbits of all major planets (Mercury through Pluto) plus the Sun over
~154 years (1899–2053) using JPL DE421 ephemeris data. The coordinate system is
rotated so Earth's orbit lies in the (x, y) plane.

```bash
uv run plot.py
```

**Output:** `solar_system.png`

---

### `plot2.py` — Sun & Moon on the Sky

Computes and plots the apparent positions of the Sun and Moon as seen from a
given observer location. Shows them as correctly-sized discs on an
altitude/azimuth chart — useful for visualizing solar eclipses.

Hardcoded locations:
- Los Alamos, NM (current time)

```bash
uv run plot2.py
```

**Output:** `e1.pdf`

---

### `plot3.py` — Solar Eclipse Simulation

Advanced version of `plot2.py` with custom light-travel-time correction,
aberration, and gravitational deflection computed from scratch (not using
Skyfield's built-in `.apparent()` method). Includes validation against
reference positions for known eclipse events.

Hardcoded eclipse events:
- Los Alamos, NM — current time
- Los Alamos, NM — Oct 14, 2023 annular eclipse
- Albuquerque, NM — Oct 14, 2023 annular eclipse
- Fredericksburg, TX — Apr 8, 2024 total eclipse

```bash
uv run plot3.py
```

**Output:** `e1.pdf`, `e2.pdf`, `e2b.pdf`, `e3.pdf`

---

### `venus_evening_chart.py` — Venus at Sunrise

Tracks Venus's position on the sky at each sunrise from October 2023 through
April 2024, as seen from Los Alamos, NM. Marker size reflects apparent
magnitude (brighter = larger). Labels show dates and month names.

```bash
uv run venus_evening_chart.py
```

**Output:** `venus_morning_chart.png`
