#!/usr/bin/env python3
"""
Compute Sun and Moon apparent positions (altitude, azimuth) from scratch.

Uses only:
  - jplephem: for reading DE421.bsp (Chebyshev polynomial interpolation)
  - numpy: linear algebra
  - IAU 2000A nutation coefficient tables (from nutation.npz)

All algorithms are implemented from first principles following:
  - USNO Circular 179
  - IERS Conventions (2010)
  - Capitaine et al. (2003)

The results must agree exactly with Skyfield's computation chain.
"""

import numpy as np
from numpy import (sin, cos, sqrt, arcsin, arctan2, arccos, fmod, dot,
                   einsum, zeros, array, where, abs, minimum, clip, outer)
import os
import sys

# SPK backend: use --pure flag or SPK_BACKEND=pure for the pure-Python reader
_use_pure = ('--pure' in sys.argv) or (os.environ.get('SPK_BACKEND') == 'pure')
if _use_pure:
    from spk_reader import SPK
else:
    from jplephem.spk import SPK

# ═══════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════

T0 = 2451545.0                      # J2000.0 epoch (Julian Date)
tau = 2 * np.pi
DAY_S = 86400.0                     # seconds per day
C_AUDAY = 173.14463267424034        # speed of light (AU/day)
AU_M = 149597870700                 # astronomical unit (meters)
AU_KM = AU_M / 1000.0
ERAD = 6378136.6                    # Earth equatorial radius (m) for deflection
ANGVEL = 7.2921150e-5               # Earth rotation rate (rad/s)
ASEC2RAD = 4.848136811095359935899141e-6  # arcsecond to radians
ASEC360 = 1296000.0                 # arcseconds in a full circle
GS = 1.32712440017987e+20           # GM_sun (m^3/s^2)
C_SI = 299792458.0                  # speed of light (m/s)

# WGS84 ellipsoid
WGS84_RADIUS = 6378137.0            # equatorial radius (m)
WGS84_INVF = 298.257223563          # inverse flattening

# Reciprocal masses (solar mass / body mass) for deflection
RMASSES = {
    'sun': 1.0, 'jupiter': 1047.3486, 'saturn': 3497.898,
    'moon': 27068700.387534, 'venus': 408523.71,
    'uranus': 22902.98, 'neptune': 19412.24,
    'earth': 332946.050895,
}

# Deflector names in priority order, and their SPK target codes
DEFLECTORS = ['sun', 'jupiter', 'saturn']
DEFLECTOR_TARGETS = {'sun': 10, 'jupiter': 5, 'saturn': 6, 'earth': 3}


# ═══════════════════════════════════════════════════════════════════════
#  Load IAU 2000A nutation coefficient tables
# ═══════════════════════════════════════════════════════════════════════

_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'skyfield', 'data')
_nut = np.load(os.path.join(_data_dir, 'nutation.npz'))

_ke0_t = _nut['ke0_t']
_ke1 = _nut['ke1']
_lunisolar_lon = _nut['lunisolar_longitude_coefficients']
_lunisolar_obl = _nut['lunisolar_obliquity_coefficients']
_nals_t = _nut['nals_t']
_napl_t = _nut['napl_t']
_nut_lon = _nut['nutation_coefficients_longitude']
_nut_obl = _nut['nutation_coefficients_obliquity']
_se0_t_0 = _nut['se0_t_0']
_se0_t_1 = _nut['se0_t_1']
_se1_0 = -0.87e-6  # hardcoded in Skyfield
_se1_1 = +0.00e-6

_TENTH_USEC_2_RAD = ASEC2RAD / 1e7

# Fundamental argument polynomial coefficients (Simon et al. 1994)
_fa0, _fa1, _fa2, _fa3, _fa4 = array((
    (485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470),
    (1287104.79305, 129596581.0481, -0.5532, 0.000136, -0.00001149),
    (335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417),
    (1072260.70369, 1602961601.2090, -6.3706, 0.006593, -0.00003169),
    (450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939),
)).T[:, :, None]

# Planetary anomaly constants and coefficients
_anomaly_constant, _anomaly_coefficient = array((
    (2.35555598, 8328.6914269554),
    (6.24006013, 628.301955),
    (1.627905234, 8433.466158131),
    (5.198466741, 7771.3771468121),
    (2.18243920, -33.757045),
    (4.402608842, 2608.7903141574),
    (3.176146697, 1021.3285546211),
    (1.753470314, 628.3075849991),
    (6.203480913, 334.0612426700),
    (0.599546497, 52.9690962641),
    (0.874016757, 21.3299104960),
    (5.481293871, 7.4781598567),
    (5.321159000, 3.8127774000),
    (0.02438175, 0.00000538691),
)).T


# ═══════════════════════════════════════════════════════════════════════
#  Rotation matrices
# ═══════════════════════════════════════════════════════════════════════

def rot_x(angle):
    c, s = cos(angle), sin(angle)
    return array([[1, 0, 0], [0, c, -s], [0, s, c]])

def rot_y(angle):
    c, s = cos(angle), sin(angle)
    return array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def rot_z(angle):
    c, s = cos(angle), sin(angle)
    return array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


# ═══════════════════════════════════════════════════════════════════════
#  Time conversions
# ═══════════════════════════════════════════════════════════════════════

def julian_day(year, month, day):
    """Compute Julian Day Number for a calendar date (integer, at noon)."""
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + A // 4
    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524

def utc_to_tt(jd_whole, utc_frac, leap_seconds=37):
    """Convert UTC Julian Date to TT Julian Date (whole + fraction)."""
    tt_frac = utc_frac + (leap_seconds + 32.184) / DAY_S
    return jd_whole, tt_frac

def tdb_minus_tt(jd_whole, tt_frac):
    """TDB - TT in seconds (USNO Circular 179, eq. 2.6)."""
    t = (jd_whole - T0 + tt_frac) / 36525.0
    return (0.001657 * sin(628.3076 * t + 6.2401)
          + 0.000022 * sin(575.3385 * t + 4.2970)
          + 0.000014 * sin(1256.6152 * t + 6.1969)
          + 0.000005 * sin(606.9777 * t + 4.0212)
          + 0.000005 * sin(52.9691 * t + 0.4444)
          + 0.000002 * sin(21.3299 * t + 5.5431)
          + 0.000010 * t * sin(628.3076 * t + 4.2490))

def tt_to_tdb(jd_whole, tt_frac):
    """Convert TT to TDB (whole + fraction)."""
    tdb_frac = tt_frac + tdb_minus_tt(jd_whole, tt_frac) / DAY_S
    return jd_whole, tdb_frac

def tt_to_ut1(jd_whole, tt_frac, delta_t):
    """Convert TT to UT1 (whole + fraction). delta_t = TT - UT1 in seconds."""
    ut1_frac = tt_frac - delta_t / DAY_S
    return jd_whole, ut1_frac


# ═══════════════════════════════════════════════════════════════════════
#  Ephemeris (SPK) queries
# ═══════════════════════════════════════════════════════════════════════

def spk_position(kernel, jd_whole, jd_frac, center, target):
    """Get position [km] from SPK kernel using two-part JD for precision."""
    return array(kernel[center, target].compute(jd_whole, jd_frac))

def spk_position_and_velocity(kernel, jd_whole, jd_frac, center, target):
    """Get position [km] and velocity [km/day] using two-part JD."""
    pv = kernel[center, target].compute_and_differentiate(jd_whole, jd_frac)
    return array(pv[0]), array(pv[1])

def earth_position_au(kernel, jd_whole, jd_frac):
    """Earth (399) position and velocity relative to SSB, in AU and AU/day."""
    emb_p, emb_v = spk_position_and_velocity(kernel, jd_whole, jd_frac, 0, 3)
    e_p, e_v = spk_position_and_velocity(kernel, jd_whole, jd_frac, 3, 399)
    return (emb_p + e_p) / AU_KM, (emb_v + e_v) / AU_KM

def sun_position_and_velocity_au(kernel, jd_whole, jd_frac):
    """Sun (10) position (AU) and velocity (AU/day) relative to SSB."""
    p, v = spk_position_and_velocity(kernel, jd_whole, jd_frac, 0, 10)
    return p / AU_KM, v / AU_KM

def moon_position_and_velocity_au(kernel, jd_whole, jd_frac):
    """Moon (301) position (AU) and velocity (AU/day) relative to SSB."""
    emb_p, emb_v = spk_position_and_velocity(kernel, jd_whole, jd_frac, 0, 3)
    m_p, m_v = spk_position_and_velocity(kernel, jd_whole, jd_frac, 3, 301)
    return (emb_p + m_p) / AU_KM, (emb_v + m_v) / AU_KM

def body_ssb_position_au(kernel, jd_whole, jd_frac, target):
    """Get barycentric position of a body (by SPK target code) in AU."""
    if target == 10:  # Sun
        return spk_position(kernel, jd_whole, jd_frac, 0, 10) / AU_KM
    elif target == 3:  # Earth
        emb = spk_position(kernel, jd_whole, jd_frac, 0, 3)
        e = spk_position(kernel, jd_whole, jd_frac, 3, 399)
        return (emb + e) / AU_KM
    else:
        return spk_position(kernel, jd_whole, jd_frac, 0, target) / AU_KM


# ═══════════════════════════════════════════════════════════════════════
#  WGS84 → ITRS
# ═══════════════════════════════════════════════════════════════════════

def wgs84_to_itrs_au(lat_deg, lon_deg, elevation_m):
    """Convert WGS84 lat/lon/elevation to ITRS position in AU."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    f = 1.0 / WGS84_INVF
    omf2 = (1.0 - f) ** 2

    sinlat = sin(lat)
    coslat = cos(lat)

    c = 1.0 / sqrt(coslat**2 + sinlat**2 * omf2)
    s = omf2 * c

    xy = (WGS84_RADIUS * c + elevation_m) * coslat
    x = xy * cos(lon)
    y = xy * sin(lon)
    z = (WGS84_RADIUS * s + elevation_m) * sinlat

    pos_m = array([x, y, z])
    return pos_m / AU_M

def itrs_velocity_au_per_day(itrs_pos_au):
    """Velocity of observer in ITRS due to Earth rotation (AU/day)."""
    x, y, z = itrs_pos_au
    return ANGVEL * DAY_S * array([-y, x, 0.0])


# ═══════════════════════════════════════════════════════════════════════
#  Precession (Capitaine et al. 2003)
# ═══════════════════════════════════════════════════════════════════════

def compute_precession(jd_tdb):
    """Compute the precession rotation matrix P for epoch jd_tdb."""
    eps0 = 84381.406  # obliquity at J2000 in arcseconds
    t = (jd_tdb - T0) / 36525.0

    psia = ((((-0.0000000951 * t
               + 0.000132851) * t
               - 0.00114045) * t
               - 1.0790069) * t
               + 5038.481507) * t

    omegaa = ((((+0.0000003337 * t
                 - 0.000000467) * t
                 - 0.00772503) * t
                 + 0.0512623) * t
                 - 0.025754) * t + eps0

    chia = ((((-0.0000000560 * t
               + 0.000170663) * t
               - 0.00121197) * t
               - 2.3814292) * t
               + 10.556403) * t

    eps0 = eps0 * ASEC2RAD
    psia = psia * ASEC2RAD
    omegaa = omegaa * ASEC2RAD
    chia = chia * ASEC2RAD

    sa = sin(eps0)
    ca = cos(eps0)
    sb = sin(-psia)
    cb = cos(-psia)
    sc = sin(-omegaa)
    cc = cos(-omegaa)
    sd = sin(chia)
    cd = cos(chia)

    return array(((cd * cb - sb * sd * cc,
                   cd * sb * ca + sd * cc * cb * ca - sa * sd * sc,
                   cd * sb * sa + sd * cc * cb * sa + ca * sd * sc),
                  (-sd * cb - sb * cd * cc,
                   -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc,
                   -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc),
                  (sb * sc,
                   -sc * cb * ca - sa * cc,
                   -sc * cb * sa + cc * ca)))


# ═══════════════════════════════════════════════════════════════════════
#  Nutation (IAU 2000A)
# ═══════════════════════════════════════════════════════════════════════

def fundamental_arguments(t):
    """5 fundamental arguments of Sun and Moon (radians). t in TDB centuries."""
    a = _fa4 * t
    a += _fa3; a *= t
    a += _fa2; a *= t
    a += _fa1; a *= t
    a += _fa0
    fmod(a, ASEC360, out=a)
    a *= ASEC2RAD
    return a[:, 0]

def iau2000a(jd_tt):
    """IAU 2000A nutation. Returns (dpsi, deps) in tenths of micro-arcseconds."""
    t = (jd_tt - T0) / 36525.0
    a = fundamental_arguments(t)

    # Luni-solar nutation
    arg = _nals_t.dot(a)
    sarg = sin(arg)
    carg = cos(arg)

    dpsi = dot(sarg, _lunisolar_lon[:, 0])
    dpsi += dot(sarg, _lunisolar_lon[:, 1]) * t
    dpsi += dot(carg, _lunisolar_lon[:, 2])

    deps = dot(carg, _lunisolar_obl[:, 0])
    deps += dot(carg, _lunisolar_obl[:, 1]) * t
    deps += dot(sarg, _lunisolar_obl[:, 2])

    # Planetary nutation
    a_plan = t * _anomaly_coefficient + _anomaly_constant
    a_plan[-1] *= t
    arg_plan = _napl_t.dot(a_plan)
    sarg_p = sin(arg_plan)
    carg_p = cos(arg_plan)

    dpsi += dot(sarg_p, _nut_lon[:, 0])
    dpsi += dot(carg_p, _nut_lon[:, 1])

    deps += dot(sarg_p, _nut_obl[:, 0])
    deps += dot(carg_p, _nut_obl[:, 1])

    return dpsi, deps

def nutation_angles_radians(jd_tt):
    """Return (d_psi, d_eps) nutation angles in radians."""
    dpsi, deps = iau2000a(jd_tt)
    return dpsi * _TENTH_USEC_2_RAD, deps * _TENTH_USEC_2_RAD

def mean_obliquity_rad(jd_tdb):
    """Mean obliquity of the ecliptic (radians). Capitaine et al. (2003)."""
    t = (jd_tdb - T0) / 36525.0
    epsilon_asec = ((((-0.0000000434 * t
                       - 0.000000576) * t
                       + 0.00200340) * t
                       - 0.0001831) * t
                       - 46.836769) * t + 84381.406
    return epsilon_asec * ASEC2RAD

def build_nutation_matrix(mean_ob, true_ob, d_psi):
    """3x3 nutation rotation matrix."""
    cobm = cos(mean_ob)
    sobm = sin(mean_ob)
    cobt = cos(true_ob)
    sobt = sin(true_ob)
    cpsi = cos(d_psi)
    spsi = sin(d_psi)

    return array(((cpsi, -spsi * cobm, -spsi * sobm),
                  (spsi * cobt, cpsi * cobm * cobt + sobm * sobt,
                   cpsi * sobm * cobt - cobm * sobt),
                  (spsi * sobt, cpsi * cobm * sobt - sobm * cobt,
                   cpsi * sobm * sobt + cobm * cobt)))


# ═══════════════════════════════════════════════════════════════════════
#  Equation of the equinoxes (complementary terms)
# ═══════════════════════════════════════════════════════════════════════

def equation_of_equinoxes_complementary_terms(jd_tt):
    """Complementary terms of the equation of the equinoxes (radians)."""
    t = (jd_tt - T0) / 36525.0
    fa = zeros(14)

    fa[0] = ((485868.249036 + (715923.2178 + (31.8792 + (0.051635 +
              (-0.00024470) * t) * t) * t) * t) * ASEC2RAD
              + (1325.0 * t % 1.0) * tau)
    fa[1] = ((1287104.793048 + (1292581.0481 + (-0.5532 + (0.000136 +
              (-0.00001149) * t) * t) * t) * t) * ASEC2RAD
              + (99.0 * t % 1.0) * tau)
    fa[2] = ((335779.526232 + (295262.8478 + (-12.7512 + (-0.001037 +
              (0.00000417) * t) * t) * t) * t) * ASEC2RAD
              + (1342.0 * t % 1.0) * tau)
    fa[3] = ((1072260.703692 + (1105601.2090 + (-6.3706 + (0.006593 +
              (-0.00003169) * t) * t) * t) * t) * ASEC2RAD
              + (1236.0 * t % 1.0) * tau)
    fa[4] = ((450160.398036 + (-482890.5431 + (7.4722 + (0.007702 +
              (-0.00005939) * t) * t) * t) * t) * ASEC2RAD
              + (-5.0 * t % 1.0) * tau)
    fa[5] = 4.402608842 + 2608.7903141574 * t
    fa[6] = 3.176146697 + 1021.3285546211 * t
    fa[7] = 1.753470314 + 628.3075849991 * t
    fa[8] = 6.203480913 + 334.0612426700 * t
    fa[9] = 0.599546497 + 52.9690962641 * t
    fa[10] = 0.874016757 + 21.3299104960 * t
    fa[11] = 5.481293872 + 7.4781598567 * t
    fa[12] = 5.311886287 + 3.8133035638 * t
    fa[13] = (0.024381750 + 0.00000538691 * t) * t

    fa %= tau

    a = _ke1.dot(fa)
    c_terms = _se1_0 * sin(a)
    c_terms += _se1_1 * cos(a)
    c_terms *= t

    a = _ke0_t.dot(fa)
    c_terms += _se0_t_0.dot(sin(a))
    c_terms += _se0_t_1.dot(cos(a))

    c_terms *= ASEC2RAD
    return c_terms


# ═══════════════════════════════════════════════════════════════════════
#  Earth rotation and sidereal time
# ═══════════════════════════════════════════════════════════════════════

def earth_rotation_angle(jd_ut1_whole, jd_ut1_frac):
    """Earth Rotation Angle (fraction of full rotation). IAU 2000 Res B1.8."""
    th = 0.7790572732640 + 0.00273781191135448 * (jd_ut1_whole - T0 + jd_ut1_frac)
    return (th % 1.0 + jd_ut1_whole % 1.0 + jd_ut1_frac) % 1.0

def greenwich_mean_sidereal_time(jd_ut1_whole, jd_ut1_frac, jd_tdb):
    """GMST in hours. Circular 179, Section 2.6.2."""
    theta = earth_rotation_angle(jd_ut1_whole, jd_ut1_frac)
    t = (jd_tdb - T0) / 36525.0

    st = (0.014506 +
          ((((-0.0000000368 * t
              - 0.000029956) * t
              - 0.00000044) * t
              + 1.3915817) * t
              + 4612.156534) * t)

    return (st / 54000.0 + theta * 24.0) % 24.0

def greenwich_apparent_sidereal_time(gmst, d_psi, mean_ob, jd_tt):
    """GAST in hours."""
    c_terms = equation_of_equinoxes_complementary_terms(jd_tt)
    eq_eq = d_psi * cos(mean_ob) + c_terms
    return (gmst + eq_eq / tau * 24.0) % 24.0


# ═══════════════════════════════════════════════════════════════════════
#  ICRS-to-J2000 bias matrix (IERS 2003 Conventions, Chapter 5)
# ═══════════════════════════════════════════════════════════════════════

def icrs_to_j2000_bias():
    """Build the frame bias rotation matrix B."""
    xi0 = -0.0166170 * ASEC2RAD
    eta0 = -0.0068192 * ASEC2RAD
    da0 = -0.01460 * ASEC2RAD

    yx = -da0
    zx = xi0
    xy = da0
    zy = eta0
    xz = -xi0
    yz = -eta0

    xx = 1.0 - 0.5 * (yx * yx + zx * zx)
    yy = 1.0 - 0.5 * (yx * yx + zy * zy)
    zz = 1.0 - 0.5 * (zy * zy + zx * zx)

    return array(((xx, xy, xz), (yx, yy, yz), (zx, zy, zz)))


# ═══════════════════════════════════════════════════════════════════════
#  Frame transformations
# ═══════════════════════════════════════════════════════════════════════

def compute_M(jd_tt, jd_tdb):
    """Compute M = N × P × B (ICRS → true equator & equinox of date)."""
    B = icrs_to_j2000_bias()
    P = compute_precession(jd_tdb)
    d_psi, d_eps = nutation_angles_radians(jd_tt)
    mean_ob = mean_obliquity_rad(jd_tdb)
    true_ob = mean_ob + d_eps
    N = build_nutation_matrix(mean_ob, true_ob, d_psi)
    M = N @ P @ B
    return M, d_psi, d_eps, mean_ob

def itrs_rotation(gast_hours, M):
    """ITRS rotation matrix: R_ITRS = Rz(-GAST) × M (no polar motion)."""
    return rot_z(-gast_hours * tau / 24.0) @ M

def altaz_rotation(lat_rad, lon_rad, R_itrs):
    """Alt-az rotation: R_altaz = R_latlon × R_ITRS.

    R_latlon = rot_y(lat)[::-1] × rot_z(-lon), following Skyfield convention
    where [::-1] reverses the row order of the latitude rotation.
    """
    R_lat = rot_y(lat_rad)[::-1]  # Skyfield convention
    R_latlon = R_lat @ rot_z(-lon_rad)
    return R_latlon @ R_itrs


# ═══════════════════════════════════════════════════════════════════════
#  Vector utilities
# ═══════════════════════════════════════════════════════════════════════

def length_of(xyz):
    return sqrt((xyz * xyz).sum(axis=0))

def dots(a, b):
    return einsum('a...,a...', a, b)

def to_spherical(xyz):
    """Convert xyz to (r, elevation, azimuth)."""
    r = length_of(xyz)
    x, y, z = xyz
    eps = np.finfo(np.float64).tiny
    theta = arcsin(z / (r + eps))       # elevation (altitude)
    phi = arctan2(y, x) % tau           # azimuth
    return r, theta, phi


# ═══════════════════════════════════════════════════════════════════════
#  Light-travel-time correction
# ═══════════════════════════════════════════════════════════════════════

def correct_for_light_travel_time(observer_pos_au, observer_vel_au_per_d,
                                  kernel, jd_tdb_whole, jd_tdb_frac,
                                  target_pos_vel_func):
    """Iteratively correct for light travel time.

    target_pos_vel_func(kernel, whole, frac) returns
    (position_au, velocity_au_per_d) relative to SSB.
    """
    t_pos, t_vel = target_pos_vel_func(kernel, jd_tdb_whole, jd_tdb_frac)

    distance = length_of(t_pos - observer_pos_au)
    light_time0 = 0.0

    for _ in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if abs(delta) < 1e-12:
            break
        frac2 = jd_tdb_frac - light_time
        t_pos, t_vel = target_pos_vel_func(kernel, jd_tdb_whole, frac2)
        distance = length_of(t_pos - observer_pos_au)
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')

    astrometric = t_pos - observer_pos_au
    astrometric_vel = t_vel - observer_vel_au_per_d
    return astrometric, astrometric_vel, light_time


# ═══════════════════════════════════════════════════════════════════════
#  Gravitational light deflection
# ═══════════════════════════════════════════════════════════════════════

def light_time_difference(position, observer_position):
    """Light-time difference between SSB and observer for a distant source."""
    dis = length_of(position)
    _AVOID = 1e-300
    u1 = position / (dis + _AVOID)
    return einsum('i,i', u1, observer_position) / C_AUDAY

def _add_deflection(position, observer, deflector, rmass):
    """Correct position in-place for one gravitating body's deflection."""
    _AVOID = 1e-300
    pq = observer + position - deflector
    pe = observer - deflector

    pmag = length_of(position)
    qmag = length_of(pq)
    emag = length_of(pe)

    phat = position / max(pmag, _AVOID)
    qhat = pq / max(qmag, _AVOID)
    ehat = pe / max(emag, _AVOID)

    pdotq = dots(phat, qhat)
    qdote = dots(qhat, ehat)
    edotp = dots(ehat, phat)

    if abs(edotp) > 0.99999999999:
        return

    fac1 = 2.0 * GS / (C_SI * C_SI * emag * AU_M * rmass)
    fac2 = 1.0 + qdote

    position += fac1 * (pdotq * ehat - edotp * qhat) / fac2 * pmag

def compute_limb_angle(position_au, observer_gcrs_au):
    """Angle of object above/below Earth's limb and nadir angle fraction."""
    earth_radius_au = ERAD / AU_M
    disobj = length_of(position_au)
    disobs = length_of(observer_gcrs_au)

    aprad = arcsin(min(earth_radius_au / disobs, 1.0))
    zdlim = np.pi - aprad

    coszd = dots(position_au, observer_gcrs_au) / (disobj * disobs)
    coszd = np.clip(coszd, -1.0, 1.0)
    zdobj = arccos(coszd)

    limb_angle = (zdlim - zdobj) * (180.0 / np.pi)
    nadir_angle = (np.pi - zdobj) / aprad
    return limb_angle, nadir_angle

def add_deflection(position, observer_bcrs, observer_gcrs, kernel,
                   jd_whole, jd_frac, count=3):
    """Apply gravitational deflection from Sun, Jupiter, Saturn (+Earth)."""
    tlt = length_of(position) / C_AUDAY

    tclose_whole = jd_whole
    tclose_frac = jd_frac
    for name in DEFLECTORS[:count]:
        target = DEFLECTOR_TARGETS[name]
        bposition = body_ssb_position_au(kernel, jd_whole, jd_frac, target)
        gpv = bposition - observer_bcrs
        dlt = light_time_difference(position, gpv)

        tclose_frac = jd_frac
        if dlt > 0.0:
            tclose_frac = jd_frac - dlt
        if tlt < dlt:
            tclose_frac = jd_frac - tlt

        bposition = body_ssb_position_au(kernel, tclose_whole, tclose_frac,
                                         target)
        _add_deflection(position, observer_bcrs, bposition, RMASSES[name])

    # Earth deflection (for observer not at geocenter)
    limb_angle, nadir_angle = compute_limb_angle(position, observer_gcrs)
    if nadir_angle >= 0.8:
        target = DEFLECTOR_TARGETS['earth']
        bposition = body_ssb_position_au(kernel, tclose_whole, tclose_frac,
                                         target)
        _add_deflection(position, observer_bcrs, bposition, RMASSES['earth'])


# ═══════════════════════════════════════════════════════════════════════
#  Aberration
# ═══════════════════════════════════════════════════════════════════════

def add_aberration(position, velocity, light_time):
    """Correct position in-place for stellar aberration."""
    _AVOID = 1e-300
    p1mag = light_time * C_AUDAY
    vemag = length_of(velocity)
    beta = vemag / C_AUDAY
    dot_pv = dots(position, velocity)

    cosd = dot_pv / (p1mag * vemag + _AVOID)
    gammai = sqrt(1.0 - beta * beta)
    p = beta * cosd
    q = (1.0 + p / (1.0 + gammai)) * light_time
    r = 1.0 + p

    position *= gammai
    position += q * velocity
    position /= r


# ═══════════════════════════════════════════════════════════════════════
#  Main computation
# ═══════════════════════════════════════════════════════════════════════

def compute_altaz(kernel, lat_deg, lon_deg, elev_m,
                  utc_year, utc_month, utc_day,
                  utc_hour, utc_minute, utc_second,
                  tz_offset_hours, delta_t, leap_seconds):
    """
    Compute Sun and Moon altitude, azimuth, and angular radius.

    Parameters:
        kernel: opened SPK file
        lat_deg, lon_deg, elev_m: observer location (WGS84)
        utc_*: local time components
        tz_offset_hours: UTC offset (e.g., -5 for CDT)
        delta_t: TT - UT1 in seconds (from IERS data)
        leap_seconds: TAI - UTC in seconds
    """
    # ── Convert local time to UTC ──
    utc_h = utc_hour - tz_offset_hours
    utc_day_offset = utc_day
    if utc_h >= 24:
        utc_h -= 24
        utc_day_offset += 1

    jd_int = julian_day(utc_year, utc_month, utc_day_offset)
    utc_frac = (utc_h * 3600 + utc_minute * 60 + utc_second) / DAY_S - 0.5

    # ── Time scales ──
    jd_whole, tt_frac = utc_to_tt(jd_int, utc_frac, leap_seconds)
    jd_tt = jd_whole + tt_frac
    _, tdb_frac = tt_to_tdb(jd_whole, tt_frac)
    jd_tdb = jd_whole + tdb_frac
    _, ut1_frac = tt_to_ut1(jd_whole, tt_frac, delta_t)

    print(f"  JD (TT):  {jd_tt:.15f}")
    print(f"  JD (TDB): {jd_tdb:.15f}")
    print(f"  delta_T:  {delta_t} s")

    # ── Precession, nutation, sidereal time ──
    M, d_psi, d_eps, mean_ob = compute_M(jd_tt, jd_tdb)
    gmst = greenwich_mean_sidereal_time(jd_whole, ut1_frac, jd_tdb)
    gast = greenwich_apparent_sidereal_time(gmst, d_psi, mean_ob, jd_tt)
    R_itrs = itrs_rotation(gast, M)

    print(f"  GMST: {gmst:.15f} h")
    print(f"  GAST: {gast:.15f} h")

    # ── Observer position ──
    itrs_pos = wgs84_to_itrs_au(lat_deg, lon_deg, elev_m)
    itrs_vel = itrs_velocity_au_per_day(itrs_pos)

    # ITRS → GCRS
    RT = R_itrs.T
    observer_gcrs = RT @ itrs_pos
    observer_vel_gcrs = RT @ itrs_vel

    # ── Earth barycentric position ──
    earth_pos, earth_vel = earth_position_au(kernel, jd_whole, tdb_frac)

    # Observer barycentric
    obs_bcrs_pos = earth_pos + observer_gcrs
    obs_bcrs_vel = earth_vel + observer_vel_gcrs

    print(f"  Observer BCRS: {obs_bcrs_pos}")
    print(f"  Observer GCRS: {observer_gcrs}")

    # ── Alt-az rotation matrix ──
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)
    R_altaz = altaz_rotation(lat_rad, lon_rad, R_itrs)

    # ── Compute for Sun ──
    # IAU 2015 solar radius (used by JPL Horizons)
    solar_radius_km = 695700.0

    sun_astro, sun_astro_vel, sun_lt = correct_for_light_travel_time(
        obs_bcrs_pos, obs_bcrs_vel, kernel, jd_whole, tdb_frac,
        sun_position_and_velocity_au)

    print(f"\n  Sun astrometric: {sun_astro}")
    print(f"  Sun light time:  {sun_lt:.15f} days")

    # Apply deflection and aberration
    sun_apparent = sun_astro.copy()
    add_deflection(sun_apparent, obs_bcrs_pos, observer_gcrs, kernel,
                   jd_whole, tdb_frac)
    add_aberration(sun_apparent, obs_bcrs_vel, sun_lt)

    # Convert to alt-az
    sun_altaz = R_altaz @ sun_apparent
    sun_dist, sun_alt, sun_az = to_spherical(sun_altaz)
    sun_ang_diam = 2.0 * np.arcsin(solar_radius_km / (sun_dist * AU_KM))
    sun_ang_diam_arcsec = np.degrees(sun_ang_diam) * 3600.0

    # ── Compute for Moon ──
    # IAU moon radius (used by JPL Horizons)
    moon_radius_km = 1737.4

    moon_astro, moon_astro_vel, moon_lt = correct_for_light_travel_time(
        obs_bcrs_pos, obs_bcrs_vel, kernel, jd_whole, tdb_frac,
        moon_position_and_velocity_au)

    print(f"  Moon astrometric: {moon_astro}")
    print(f"  Moon light time:  {moon_lt:.15e} days")

    moon_apparent = moon_astro.copy()
    add_deflection(moon_apparent, obs_bcrs_pos, observer_gcrs, kernel,
                   jd_whole, tdb_frac)
    add_aberration(moon_apparent, obs_bcrs_vel, moon_lt)

    moon_altaz = R_altaz @ moon_apparent
    moon_dist, moon_alt, moon_az = to_spherical(moon_altaz)
    moon_ang_diam = 2.0 * np.arcsin(moon_radius_km / (moon_dist * AU_KM))
    moon_ang_diam_arcsec = np.degrees(moon_ang_diam) * 3600.0

    # ── Print results ──
    sun_alt_deg = np.degrees(sun_alt)
    sun_az_deg = np.degrees(sun_az)
    moon_alt_deg = np.degrees(moon_alt)
    moon_az_deg = np.degrees(moon_az)

    print(f"\nSun:")
    print(f"  Altitude:  {sun_alt_deg:.6f}°")
    print(f"  Azimuth:   {sun_az_deg:.6f}°")
    print(f"  Ang-diam:  {sun_ang_diam_arcsec:.3f}\"")
    print(f"  delta:     {sun_dist:.14f} AU")
    print(f"Moon:")
    print(f"  Altitude:  {moon_alt_deg:.6f}°")
    print(f"  Azimuth:   {moon_az_deg:.6f}°")
    print(f"  Ang-diam:  {moon_ang_diam_arcsec:.3f}\"")
    print(f"  delta:     {moon_dist:.14f} AU")

    return (sun_alt_deg, sun_az_deg, sun_ang_diam_arcsec, sun_dist,
            moon_alt_deg, moon_az_deg, moon_ang_diam_arcsec, moon_dist)


# ═══════════════════════════════════════════════════════════════════════
#  Entry point
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    backend = 'spk_reader (pure Python)' if _use_pure else 'jplephem'
    bsp_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'de440s.bsp')
    kernel = SPK.open(bsp_path)

    print("=" * 65)
    print("40°N, Greenwich — 2025 January 1, 12:00 UTC")
    print(f"SPK backend: {backend}")
    print("=" * 65)

    result = compute_altaz(
        kernel,
        lat_deg=40.0, lon_deg=0.0, elev_m=0,
        utc_year=2025, utc_month=1, utc_day=1,
        utc_hour=12, utc_minute=0, utc_second=0,
        tz_offset_hours=0,    # UTC
        delta_t=69.14980035,  # TT - UT1 (seconds), from IERS data
        leap_seconds=37,      # TAI - UTC for 2025
    )

    # ── Validate against JPL Horizons (DE441) reference values ──
    # Horizons uses DE441, we use DE440s — small differences expected in
    # position due to different ephemeris and EOP data. Angular diameters
    # should match closely since they depend mainly on distance and radii.
    horizons = {
        'Sun alt':       (result[0], 27.036034),
        'Sun az':        (result[1], 179.049603),
        'Sun Ang-diam':  (result[2], 1950.991),
        'Sun delta':     (result[3], 0.98332708143732),
        'Moon alt':      (result[4], 21.518703),
        'Moon az':       (result[5], 157.820214),
        'Moon Ang-diam': (result[6], 1897.634),
        'Moon delta':    (result[7], 0.00252475127904),
    }

    print("\n" + "=" * 65)
    print("Comparison with JPL Horizons (DE441)")
    print("=" * 65)
    for label, (computed, expected) in horizons.items():
        diff = computed - expected
        adiff = abs(diff)
        if 'diam' in label:
            unit = '"'
            print(f"  {label:16s}: {computed:12.3f}  Horizons: {expected:12.3f}  diff: {diff:+.3f}{unit}")
        elif 'delta' in label:
            unit = ' AU'
            print(f"  {label:16s}: {computed:.14f}  Horizons: {expected:.14f}  diff: {diff:+.2e}{unit}")
        else:
            unit = '°'
            print(f"  {label:16s}: {computed:12.6f}  Horizons: {expected:12.6f}  diff: {diff:+.6f}{unit} ({adiff*3600:.3f}\")")

    print()
    print("  Note: Horizons uses DE441 ephemeris + its own EOP data.")
    print("  We use DE440s + a single delta_T value. Small alt/az")
    print("  differences (~0.1\") are expected from these sources.")

    kernel.close()
