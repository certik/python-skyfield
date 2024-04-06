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
from skyfield.positionlib import Apparent, Astrometric, build_position
from skyfield.functions import dots
from skyfield.functions import (
    _T, _to_array, _to_spherical_and_rates, angle_between, from_spherical,
    length_of, mxm, mxv, rot_z,
)
from skyfield.constants import (AU_M, ANGVEL, DAY_S, DEG2RAD, ERAD,
                        IERS_2010_INVERSE_EARTH_FLATTENING, RAD2DEG, T0, tau,
                        C_AUDAY)
from skyfield.units import Angle, AngleRate, Distance, Velocity, _interpret_angle
from skyfield.timelib import Time



# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

earth, sun, moon, venus = eph['earth'], eph['sun'], eph['moon'], eph['venus']
ts = load.timescale()

solar_radius_km = 696340.0
moon_radius_km = 1737.1
earth_radius_au = ERAD / AU_M
AU_KM = 149597870.700

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

def compute_apparent(bcrs_position, bcrs_velocity, observer_gcrs_au, t, position, ephemeris, light_time):
    target_au = position.copy()

    limb_angle, nadir_angle = compute_limb_angle(
        target_au, observer_gcrs_au)
    include_earth_deflection = nadir_angle >= 0.8

    add_deflection(target_au, bcrs_position,
                   ephemeris, t, include_earth_deflection)

    add_aberration(target_au, bcrs_velocity, light_time)

    return target_au

def length_of(xyz):
    """Given a 3-element array |xyz|, return its length.

    The three elements can be simple scalars, or the array can be two
    dimensions and offer three whole series of x, y, and z coordinates.

    """
    return sqrt((xyz * xyz).sum(axis=0))

def to_spherical(xyz):
    """Convert |xyz| to spherical coordinates (r,theta,phi).

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    Note that ``theta`` is an elevation angle measured up and down from
    the xy-plane, not a polar angle measured from the z-axis, to match
    the convention for both latitude and declination.

    """
    r = length_of(xyz)
    x, y, z = xyz
    eps = np.finfo(np.float64).tiny
    theta = arcsin(z / (r + eps))
    phi = arctan2(y, x) % tau
    return r, theta, phi

def altaz(position, R):
    """Compute (alt, az, distance) relative to the observer's horizon.
    
    `position` is in AU
    """
    position_au = np.dot(R, position)
    r_au, alt, az = to_spherical(position_au)
    return alt, az, r_au

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

def correct_for_light_travel_time(observer, target):
    """Return a light-time corrected astrometric position and velocity.

    Given an `observer` that is a `Barycentric` position somewhere in
    the solar system, compute where in the sky they will see the body
    `target`, by computing the light-time between them and figuring out
    where `target` was back when the light was leaving it that is now
    reaching the eyes or instruments of the `observer`.

    """
    t = observer.t
    ts = t.ts
    whole = t.whole
    tdb_fraction = t.tdb_fraction

    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d

    tposition, tvelocity, gcrs_position, message = target._at(t)

    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if np.max(abs(delta), initial=0.0) < 1e-12:
            break

        # We assume a light travel time of at most a couple of days.  A
        # longer light travel time would best be split into a whole and
        # fraction, for adding to the whole and fraction of TDB.
        t2 = ts.tdb_jd(whole, tdb_fraction - light_time)

        tposition, tvelocity, gcrs_position, message = target._at(t2)
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')
    return tposition - cposition, tvelocity - cvelocity, t, light_time

def observe(self, body):
    """Compute the `Astrometric` position of a body from this location.

    To compute the body's astrometric position, it is first asked
    for its position at the time `t` of this position itself.  The
    distance to the body is then divided by the speed of light to
    find how long it takes its light to arrive.  Finally, the light
    travel time is subtracted from `t` and the body is asked for a
    series of increasingly exact positions to learn where it was
    when it emitted the light that is now reaching this position.

    >>> earth.at(t).observe(mars)
    <Astrometric ICRS position and velocity at date t center=399 target=499>

    """
    p, v, t, light_time = correct_for_light_travel_time(self, body)
    astrometric = Astrometric(p, v, t, self.target, body.target)
    astrometric._ephemeris = self._ephemeris
    astrometric.center_barycentric = self
    astrometric.light_time = light_time
    return astrometric

def at(self, t):
    """At time ``t``, compute the target's position relative to the center.

    If ``t`` is an array of times, then the returned position object
    will specify as many positions as there were times.  The kind of
    position returned depends on the value of the ``center``
    attribute:

    * Solar System Barycenter: :class:`~skyfield.positionlib.Barycentric`
    * Center of the Earth: :class:`~skyfield.positionlib.Geocentric`
    * Anything else: :class:`~skyfield.positionlib.ICRF`

    """
    if not isinstance(t, Time):
        raise ValueError('please provide the at() method with a Time'
                            ' instance as its argument, instead of the'
                            ' value {0!r}'.format(t))
    p, v, gcrs_position, message = self._at(t)
    center = self.center
    position = build_position(p, v, t, center, self.target)
    position._ephemeris = self.ephemeris
    position._observer_gcrs_au = gcrs_position
    position.message = message
    return position

def compute(observer, zone, time, loc, filename,
        ref_pos=None):
    time = zone.localize(time)
    t = ts.from_datetime(time)
    print("Location:", loc)
    print(observer)
    print("Time:", t.astimezone(zone))

    at_result = at(earth + observer, t)
    obs = observe(at_result, sun)
    #apparent = obs.apparent()
    R = obs.center_barycentric._altaz_rotation
    position = compute_apparent(
        obs.center_barycentric.position.au,
        obs.center_barycentric.velocity.au_per_d,
        obs.center_barycentric._observer_gcrs_au,
        obs.t,
        obs.position.au, obs._ephemeris, obs.light_time)
    alt, az, distance = altaz(position, R)
    sun_r = asin(solar_radius_km/(distance*AU_KM))
    sun_alt = rad_to_deg(alt)
    sun_az = rad_to_deg(az)
    print(sun_alt, sun_az, sun_r)
    print("Sun:")
    print("Altitude (0-90):", rad_to_str(alt))
    print("Azimuth (0-360): ", rad_to_str(az))
    print("Radius (deg):", rad_to_str(sun_r))
    print()
    obs = (earth + observer).at(t).observe(moon)
    R = obs.center_barycentric._altaz_rotation
    position = compute_apparent(
        obs.center_barycentric.position.au,
        obs.center_barycentric.velocity.au_per_d,
        obs.center_barycentric._observer_gcrs_au,
        obs.t,
        obs.position.au, obs._ephemeris, obs.light_time)
    alt, az, distance = altaz(position, R)
    moon_r = asin(moon_radius_km/(distance*AU_KM))
    moon_alt = rad_to_deg(alt)
    moon_az = rad_to_deg(az)
    print(moon_alt, moon_az, moon_r)
    print("Moon:")
    print("Altitude (0-90):", rad_to_str(alt))
    print("Azimuth (0-360): ", rad_to_str(az))
    print("Radius (deg):", rad_to_str(moon_r))
    print()
    print()
    print()
    if ref_pos:
        eps = 1e-16
        assert abs(sun_alt-ref_pos[0]) < eps
        assert abs(sun_az-ref_pos[1]) < eps
        assert abs(sun_r-ref_pos[2]) < eps
        assert abs(moon_alt-ref_pos[3]) < eps
        assert abs(moon_az-ref_pos[4]) < eps
        assert abs(moon_r-ref_pos[5]) < eps


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
        "Los Alamos, NM", "e2.pdf",
        [35.62384236539088, 137.39322269477546, 0.004666435109108667,
         35.62491478656301, 137.40792969645054, 0.004415235323045712])


# https://www.timeanddate.com/eclipse/in/@35.16789477316579,-106.58454895019533?iso=20231014
compute(wgs84.latlon(35.16789477316579,-106.58454895019533, 1600),
        timezone('US/Mountain'), datetime(2023, 10, 14, 10, 37, 1),
        "Albuquerque, NM", "e2b.pdf",
        [36.13290440870471, 136.97575208131278, 0.00466643734720691,
         36.132703895971886, 136.9757657675733, 0.004415759291846439],
        )

# https://www.timeanddate.com/eclipse/in/@30.276139106760382,-98.872489929199233?iso=20240408
compute(wgs84.latlon(30.274167, -98.871944, 516),
        timezone('US/Central'), datetime(2024, 4, 8, 13, 35, 10),
        "Fredericksburg, TX", "e3.pdf",
        [67.31640064112162, 178.74688214701396, 0.0046479241546984506,
         67.31763737226208, 178.74983743582982, 0.004908042955826713])