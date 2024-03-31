from datetime import datetime
from math import atan, pi
import numpy as np
from pytz import timezone
from matplotlib.pylab import (plot, savefig, legend, grid, gca, scatter,
        figure, xlim, ylim)
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from skyfield.api import load, wgs84
from skyfield.units import Angle

# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

earth, sun, moon, venus = eph['earth'], eph['sun'], eph['moon'], eph['venus']
ts = load.timescale()

solar_radius_km = 696340.0
moon_radius_km = 1737.1

def compute(observer, zone, time, loc, filename):
    time = zone.localize(time)
    t = ts.from_datetime(time)
    print("Location:", loc)
    print(observer)
    print("Time:", t.astimezone(zone))

    apparent = (earth + observer).at(t).observe(sun).apparent()
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


    figure()
    plot([sun_az], [sun_alt], "oy")
    plot([moon_az], [moon_alt], ".k")
    gca().add_patch(Circle((sun_az, sun_alt), sun_r,
        color="y", alpha=0.5))
    gca().add_patch(Circle((sun_az, sun_alt), sun_r,
        ec="orange", fc="none", lw=2))
    gca().add_patch(Circle((moon_az, moon_alt), moon_r,
        ec="k", fc="none", lw=2))
    gca().set_aspect("equal")
    grid()
    xlim([sun_az-1, sun_az+1])
    ylim([sun_alt-1, sun_alt+1])

    savefig(filename)


compute(wgs84.latlon(35.90369476314685, -106.30851268194805, 2230),
        timezone('US/Mountain'), datetime.now(),
        "Los Alamos, NM", "e1.png")

# https://eclipse2024.org/2023eclipse/eclipse-cities/city/28230.html
compute(wgs84.latlon(35.90369476314685, -106.30851268194805, 2230),
        timezone('US/Mountain'), datetime(2023, 10, 14, 10, 36, 8),
        "Los Alamos, NM", "e2.png")

# https://www.timeanddate.com/eclipse/in/@35.16789477316579,-106.58454895019533?iso=20231014
compute(wgs84.latlon(35.16789477316579,-106.58454895019533, 1600),
        timezone('US/Mountain'), datetime(2023, 10, 14, 10, 37, 1),
        "Albuquerque, NM", "e2b.png")

# https://www.timeanddate.com/eclipse/in/@30.276139106760382,-98.872489929199233?iso=20240408
compute(wgs84.latlon(30.274167, -98.871944, 516),
        timezone('US/Central'), datetime(2024, 4, 8, 13, 35, 10),
        "Fredericksburg, TX", "e3.png")
