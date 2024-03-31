from datetime import datetime
from math import atan, pi
import numpy as np
from pytz import timezone
from matplotlib.pylab import (plot, savefig, legend, grid, gca, scatter,
        figure, xlim, ylim, title, xlabel, ylabel)
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from skyfield.api import load, wgs84
from skyfield.units import Angle
from skyfield.earthlib import compute_limb_angle
from skyfield.relativity import add_aberration, add_deflection
from skyfield.positionlib import Apparent



# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

earth, sun, moon, venus = eph['earth'], eph['sun'], eph['moon'], eph['venus']
ts = load.timescale()

solar_radius_km = 696340.0
moon_radius_km = 1737.1

def compute_apparent(self):
    t = self.t
    target_au = self.position.au.copy()

    cb = self.center_barycentric
    bcrs_position = cb.position.au
    bcrs_velocity = cb.velocity.au_per_d
    observer_gcrs_au = cb._observer_gcrs_au

    if observer_gcrs_au is None:
        include_earth_deflection = array((False,))
    else:
        limb_angle, nadir_angle = compute_limb_angle(
            target_au, observer_gcrs_au)
        include_earth_deflection = nadir_angle >= 0.8

    add_deflection(target_au, bcrs_position,
                   self._ephemeris, t, include_earth_deflection)

    add_aberration(target_au, bcrs_velocity, self.light_time)

    v = self.velocity.au_per_d
    if v is not None:
        pass  # TODO: how to apply aberration and deflection to velocity?

    apparent = Apparent(target_au, v, t, self.center, self.target)
    apparent.center_barycentric = self.center_barycentric
    apparent._observer_gcrs_au = observer_gcrs_au
    return apparent


def compute(observer, zone, time, loc, filename):
    time = zone.localize(time)
    t = ts.from_datetime(time)
    print("Location:", loc)
    print(observer)
    print("Time:", t.astimezone(zone))

    obs = (earth + observer).at(t).observe(sun)
    #apparent = obs.apparent()
    apparent = compute_apparent(obs)
    alt, az, distance = apparent.altaz()
    radius_angle = Angle(radians=atan(solar_radius_km/distance.km))
    sun_alt = alt.degrees
    sun_az = az.degrees
    sun_r = radius_angle.degrees
    print("Sun:")
    print("Altitude (0-90):", alt)
    print("Azimuth (0-360): ", az)
    print("Radius (deg):", radius_angle)
    print()
    apparent = (earth + observer).at(t).observe(moon).apparent()
    alt, az, distance = apparent.altaz()
    radius_angle = Angle(radians=atan(moon_radius_km/distance.km))
    moon_alt = alt.degrees
    moon_az = az.degrees
    moon_r = radius_angle.degrees
    print("Moon:")
    print("Altitude (0-90):", alt)
    print("Azimuth (0-360): ", az)
    print("Radius (deg):", radius_angle)
    print()
    print()
    print()


    figure(figsize=(5,5))
    #plot([sun_az], [sun_alt], "oy")
    #plot([moon_az], [moon_alt], ".k")
    ax = gca()
    ax.add_patch(Circle((sun_az, sun_alt), sun_r,
        ec="none", fc="orange", lw=0))
    ax.add_patch(Circle((moon_az, moon_alt), moon_r,
        ec="none", fc="k", lw=0))
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    grid()
    xlim([sun_az-1, sun_az+1])
    ylim([sun_alt-1, sun_alt+1])
    xlabel(r"Azimuth [deg]")
    ylabel(r"Altitude [deg]")
    title("Sun and Moon on the sky\n%s\n%s" % (loc, t.astimezone(zone)))
    savefig(filename)


compute(wgs84.latlon(35.90369476314685, -106.30851268194805, 2230),
        timezone('US/Mountain'), datetime.now(),
        "Los Alamos, NM", "e1.pdf")

# https://eclipse2024.org/2023eclipse/eclipse-cities/city/28230.html
compute(wgs84.latlon(35.90369476314685, -106.30851268194805, 2230),
        timezone('US/Mountain'), datetime(2023, 10, 14, 10, 36, 8),
        "Los Alamos, NM", "e2.pdf")

# https://www.timeanddate.com/eclipse/in/@35.16789477316579,-106.58454895019533?iso=20231014
compute(wgs84.latlon(35.16789477316579,-106.58454895019533, 1600),
        timezone('US/Mountain'), datetime(2023, 10, 14, 10, 37, 1),
        "Albuquerque, NM", "e2b.pdf")

# https://www.timeanddate.com/eclipse/in/@30.276139106760382,-98.872489929199233?iso=20240408
compute(wgs84.latlon(30.274167, -98.871944, 516),
        timezone('US/Central'), datetime(2024, 4, 8, 13, 35, 10),
        "Fredericksburg, TX", "e3.pdf")
