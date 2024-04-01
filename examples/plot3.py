from datetime import datetime
from math import atan, pi, asin
import numpy as np
from numpy import (abs, arcsin, arccos, arctan2, array, clip, cos,
                   minimum, pi, sin, sqrt, tan, where, zeros_like)
from pytz import timezone
from matplotlib.pylab import (plot, savefig, legend, grid, gca, scatter,
        figure, xlim, ylim, title, xlabel, ylabel)
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from skyfield.api import load, wgs84
from skyfield.units import Angle
from skyfield.relativity import add_aberration, add_deflection
from skyfield.positionlib import Apparent
from skyfield.functions import dots
from skyfield.functions import (
    _T, _to_array, _to_spherical_and_rates, angle_between, from_spherical,
    length_of, mxm, mxv, rot_z, to_spherical,
)
from skyfield.constants import (AU_M, ANGVEL, DAY_S, DEG2RAD, ERAD,
                        IERS_2010_INVERSE_EARTH_FLATTENING, RAD2DEG, T0, tau)
from skyfield.units import Angle, AngleRate, Distance, Velocity, _interpret_angle



# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

earth, sun, moon, venus = eph['earth'], eph['sun'], eph['moon'], eph['venus']
ts = load.timescale()

solar_radius_km = 696340.0
moon_radius_km = 1737.1
earth_radius_au = ERAD / AU_M

def compute_limb_angle(position_au, observer_au):
    """Determine the angle of an object above or below the Earth's limb.

    Given an object's GCRS `position_au` |xyz| vector and the position
    of an `observer_au` as a vector in the same coordinate system,
    return a tuple that provides `(limb_ang, nadir_ang)`:

    limb_angle
        Angle of observed object above (+) or below (-) limb in degrees.
    nadir_angle
        Nadir angle of observed object as a fraction of apparent radius
        of limb: <1.0 means below the limb, =1.0 means on the limb, and
        >1.0 means above the limb.

    """
    # Compute the distance to the object and the distance to the observer.

    disobj = sqrt(dots(position_au, position_au))
    disobs = sqrt(dots(observer_au, observer_au))

    # Compute apparent angular radius of Earth's limb.

    aprad = arcsin(minimum(earth_radius_au / disobs, 1.0))

    # Compute zenith distance of Earth's limb.

    zdlim = pi - aprad

    # Compute zenith distance of observed object.

    coszd = dots(position_au, observer_au) / (disobj * disobs)
    coszd = clip(coszd, -1.0, 1.0)
    zdobj = arccos(coszd)

    # Angle of object wrt limb is difference in zenith distances.

    limb_angle = (zdlim - zdobj) * RAD2DEG

    # Nadir angle of object as a fraction of angular radius of limb.

    nadir_angle = (pi - zdobj) / aprad

    return limb_angle, nadir_angle

def compute_apparent(self):
    t = self.t
    target_au = self.position.au.copy()

    cb = self.center_barycentric
    bcrs_position = cb.position.au
    bcrs_velocity = cb.velocity.au_per_d
    observer_gcrs_au = cb._observer_gcrs_au

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

def altaz(position, R):
    """Compute (alt, az, distance) relative to the observer's horizon.
    
    `position` is in AU
    """
    position_au = np.dot(R, position)
    r_au, alt, az = to_spherical(position_au)
    return alt, az, Distance(r_au)

def deg_to_int(value, places=0):
    """Decompose `value` into units, minutes, seconds, and second fractions.

    This routine prepares a value for sexagesimal display, with its
    seconds fraction expressed as an integer with `places` digits.  The
    result is a tuple of five integers:

    ``(sign [either +1 or -1], units, minutes, seconds, second_fractions)``

    The integers are properly rounded per astronomical convention so
    that, for example, given ``places=3`` the result tuple ``(1, 11, 22,
    33, 444)`` means that the input was closer to 11u 22' 33.444" than
    to either 33.443" or 33.445" in its value.

    """
    power = 10 ** places
    n = int((power * 3600 * value + 0.5) // 1.0)
    sign = np.sign(n)
    n, fraction = divmod(abs(n), power)
    n, seconds = divmod(n, 60)
    n, minutes = divmod(n, 60)
    return sign, n, minutes, seconds, fraction

def rad_to_deg(value):
    return value*180/pi

def deg_to_str(value):
    sign, d, m, s, sf = deg_to_int(value, 1)
    if sign >= 0:
        st = ""
    else:
        st = "-"

    return f"{st}{d:02}° {m:02}' {s:02}.{sf}\""

def rad_to_str(value):
    return deg_to_str(rad_to_deg(value))

def compute(observer, zone, time, loc, filename):
    time = zone.localize(time)
    t = ts.from_datetime(time)
    print("Location:", loc)
    print(observer)
    print("Time:", t.astimezone(zone))

    obs = (earth + observer).at(t).observe(sun)
    #apparent = obs.apparent()
    apparent = compute_apparent(obs)
    alt, az, distance = altaz(apparent.position.au, apparent.center_barycentric._altaz_rotation)
    sun_r = asin(solar_radius_km/distance.km)
    sun_alt = rad_to_deg(alt)
    sun_az = rad_to_deg(az)
    print("Sun:")
    print("Altitude (0-90):", rad_to_str(alt))
    print("Azimuth (0-360): ", rad_to_str(az))
    print("Radius (deg):", rad_to_str(sun_r))
    print()
    apparent = (earth + observer).at(t).observe(moon).apparent()
    alt, az, distance = altaz(apparent.position.au, apparent.center_barycentric._altaz_rotation)
    moon_r = asin(moon_radius_km/distance.km)
    moon_alt = rad_to_deg(alt)
    moon_az = rad_to_deg(az)
    print("Moon:")
    print("Altitude (0-90):", rad_to_str(alt))
    print("Azimuth (0-360): ", rad_to_str(az))
    print("Radius (deg):", rad_to_str(moon_r))
    print()
    print()
    print()


    figure(figsize=(5,5))
    #plot([sun_az], [sun_alt], "oy")
    #plot([moon_az], [moon_alt], ".k")
    ax = gca()
    ax.add_patch(Circle((sun_az, sun_alt), rad_to_deg(sun_r),
        ec="none", fc="orange", lw=0))
    ax.add_patch(Circle((moon_az, moon_alt), rad_to_deg(moon_r),
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
