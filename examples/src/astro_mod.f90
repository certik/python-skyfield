! ─────────────────────────────────────────────────────────────────────
!  astro_mod — astronomical computations (time, precession, nutation,
!              sidereal time, observer position, light travel, deflection,
!              aberration, alt/az).  Extracted from generate_observations.f90.
! ─────────────────────────────────────────────────────────────────────
module astro_mod
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  implicit none
contains

  ! ── Julian Day number ──
  function julian_day(year, month, day) result(jd)
    integer, intent(in) :: year, month, day
    integer :: jd, y, m, A, B
    y = year; m = month
    if (m <= 2) then
      y = y - 1
      m = m + 12
    end if
    A = y / 100
    B = 2 - A + A / 4
    jd = int(365.25_dp * (y + 4716)) + int(30.6001_dp * (m + 1)) + day + B - 1524
  end function

  ! ── UTC → TT ──
  subroutine utc_to_tt(jd_whole, utc_frac, leap_sec, tt_frac_out)
    real(dp), intent(in)  :: jd_whole, utc_frac
    integer, intent(in)   :: leap_sec
    real(dp), intent(out) :: tt_frac_out
    tt_frac_out = utc_frac + (real(leap_sec, dp) + 32.184_dp) / DAY_S
  end subroutine

  ! ── TDB − TT (USNO Circular 179, eq 2.6) ──
  function tdb_minus_tt(jd_whole, tt_frac) result(dt)
    real(dp), intent(in) :: jd_whole, tt_frac
    real(dp) :: dt, t
    t = (jd_whole - T0 + tt_frac) / 36525.0_dp
    dt = 0.001657_dp * sin(628.3076_dp * t + 6.2401_dp) &
       + 0.000022_dp * sin(575.3385_dp * t + 4.2970_dp) &
       + 0.000014_dp * sin(1256.6152_dp * t + 6.1969_dp) &
       + 0.000005_dp * sin(606.9777_dp * t + 4.0212_dp) &
       + 0.000005_dp * sin(52.9691_dp * t + 0.4444_dp) &
       + 0.000002_dp * sin(21.3299_dp * t + 5.5431_dp) &
       + 0.000010_dp * t * sin(628.3076_dp * t + 4.2490_dp)
  end function

  ! ── TT → TDB ──
  subroutine tt_to_tdb(jd_whole, tt_frac, tdb_frac_out)
    real(dp), intent(in)  :: jd_whole, tt_frac
    real(dp), intent(out) :: tdb_frac_out
    tdb_frac_out = tt_frac + tdb_minus_tt(jd_whole, tt_frac) / DAY_S
  end subroutine

  ! ── TT → UT1 ──
  subroutine tt_to_ut1(jd_whole, tt_frac, delta_t, ut1_frac_out)
    real(dp), intent(in)  :: jd_whole, tt_frac, delta_t
    real(dp), intent(out) :: ut1_frac_out
    ut1_frac_out = tt_frac - delta_t / DAY_S
  end subroutine

  ! ── Precession (Capitaine 2003) ──
  function compute_precession(jd_tdb) result(P)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: P(3,3), t
    real(dp) :: eps0_as, psia, omegaa, chia
    real(dp) :: eps0, sa, ca, sb, cb, sc, cc_v, sd, cd

    t = (jd_tdb - T0) / 36525.0_dp
    eps0_as = 84381.406_dp

    psia = ((((-0.0000000951_dp * t &
               + 0.000132851_dp) * t &
               - 0.00114045_dp) * t &
               - 1.0790069_dp) * t &
               + 5038.481507_dp) * t

    omegaa = ((((+0.0000003337_dp * t &
                 - 0.000000467_dp) * t &
                 - 0.00772503_dp) * t &
                 + 0.0512623_dp) * t &
                 - 0.025754_dp) * t + eps0_as

    chia = ((((-0.0000000560_dp * t &
               + 0.000170663_dp) * t &
               - 0.00121197_dp) * t &
               - 2.3814292_dp) * t &
               + 10.556403_dp) * t

    eps0   = eps0_as * ASEC2RAD
    psia   = psia * ASEC2RAD
    omegaa = omegaa * ASEC2RAD
    chia   = chia * ASEC2RAD

    sa = sin(eps0);     ca = cos(eps0)
    sb = sin(-psia);    cb = cos(-psia)
    sc = sin(-omegaa);  cc_v = cos(-omegaa)
    sd = sin(chia);     cd = cos(chia)

    P(1,1) = cd * cb - sb * sd * cc_v
    P(1,2) = cd * sb * ca + sd * cc_v * cb * ca - sa * sd * sc
    P(1,3) = cd * sb * sa + sd * cc_v * cb * sa + ca * sd * sc
    P(2,1) = -sd * cb - sb * cd * cc_v
    P(2,2) = -sd * sb * ca + cd * cc_v * cb * ca - sa * cd * sc
    P(2,3) = -sd * sb * sa + cd * cc_v * cb * sa + ca * cd * sc
    P(3,1) = sb * sc
    P(3,2) = -sc * cb * ca - sa * cc_v
    P(3,3) = -sc * cb * sa + cc_v * ca
  end function

  ! ── Nutation matrix ──
  function build_nutation_matrix(mean_ob, true_ob, d_psi) result(N)
    real(dp), intent(in) :: mean_ob, true_ob, d_psi
    real(dp) :: N(3,3)
    real(dp) :: cobm, sobm, cobt, sobt, cpsi, spsi

    cobm = cos(mean_ob);  sobm = sin(mean_ob)
    cobt = cos(true_ob);  sobt = sin(true_ob)
    cpsi = cos(d_psi);    spsi = sin(d_psi)

    N(1,1) = cpsi
    N(1,2) = -spsi * cobm
    N(1,3) = -spsi * sobm
    N(2,1) = spsi * cobt
    N(2,2) = cpsi * cobm * cobt + sobm * sobt
    N(2,3) = cpsi * sobm * cobt - cobm * sobt
    N(3,1) = spsi * sobt
    N(3,2) = cpsi * cobm * sobt - sobm * cobt
    N(3,3) = cpsi * sobm * sobt + cobm * cobt
  end function

  ! ── ICRS → J2000 bias ──
  function icrs_to_j2000_bias() result(B)
    real(dp) :: B(3,3)
    real(dp) :: xi0, eta0, da0
    real(dp) :: yx, zx, xy, zy, xz, yz

    xi0  = -0.0166170_dp * ASEC2RAD
    eta0 = -0.0068192_dp * ASEC2RAD
    da0  = -0.01460_dp * ASEC2RAD

    yx = -da0;   zx = xi0
    xy =  da0;   zy = eta0
    xz = -xi0;   yz = -eta0

    B(1,1) = 1.0_dp - 0.5_dp * (yx*yx + zx*zx)
    B(1,2) = xy
    B(1,3) = xz
    B(2,1) = yx
    B(2,2) = 1.0_dp - 0.5_dp * (yx*yx + zy*zy)
    B(2,3) = yz
    B(3,1) = zx
    B(3,2) = zy
    B(3,3) = 1.0_dp - 0.5_dp * (zy*zy + zx*zx)
  end function

  ! ── Compute M = N × P × B ──
  subroutine compute_M(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
    real(dp), intent(in)  :: jd_tt, jd_tdb
    real(dp), intent(out) :: M(3,3), d_psi, d_eps, mean_ob
    real(dp) :: B(3,3), P(3,3), Nmat(3,3), true_ob
    real(dp) :: dpsi_raw, deps_raw

    B = icrs_to_j2000_bias()
    P = compute_precession(jd_tdb)
    call iau2000a(jd_tt, dpsi_raw, deps_raw)
    d_psi = dpsi_raw * TENTH_USEC_2_RAD
    d_eps = deps_raw * TENTH_USEC_2_RAD
    mean_ob = mean_obliquity_rad(jd_tdb)
    true_ob = mean_ob + d_eps
    Nmat = build_nutation_matrix(mean_ob, true_ob, d_psi)
    M = mat33_mul(Nmat, mat33_mul(P, B))
  end subroutine

  ! ── Earth Rotation Angle ──
  function earth_rotation_angle(jd_whole, jd_frac) result(theta)
    real(dp), intent(in) :: jd_whole, jd_frac
    real(dp) :: theta, th
    th = 0.7790572732640_dp + 0.00273781191135448_dp * &
         (jd_whole - T0 + jd_frac)
    theta = mod(mod(th, 1.0_dp) + mod(jd_whole, 1.0_dp) + jd_frac, 1.0_dp)
  end function

  ! ── GMST ──
  function greenwich_mean_sidereal_time(jd_whole, jd_frac, jd_tdb) result(gmst)
    real(dp), intent(in) :: jd_whole, jd_frac, jd_tdb
    real(dp) :: gmst, theta, t, st
    theta = earth_rotation_angle(jd_whole, jd_frac)
    t = (jd_tdb - T0) / 36525.0_dp
    st = 0.014506_dp + &
         ((((-0.0000000368_dp * t &
             - 0.000029956_dp) * t &
             - 0.00000044_dp) * t &
             + 1.3915817_dp) * t &
             + 4612.156534_dp) * t
    gmst = mod(st / 54000.0_dp + theta * 24.0_dp, 24.0_dp)
  end function

  ! ── GAST ──
  function greenwich_apparent_sidereal_time(gmst, d_psi, mean_ob, jd_tt) result(gast)
    real(dp), intent(in) :: gmst, d_psi, mean_ob, jd_tt
    real(dp) :: gast, c_terms, eq_eq
    c_terms = eq_equinox_complement(jd_tt)
    eq_eq = d_psi * cos(mean_ob) + c_terms
    gast = mod(gmst + eq_eq / TAU * 24.0_dp, 24.0_dp)
  end function

  ! ── ITRS rotation: Rz(-GAST) × M ──
  function itrs_rotation(gast_hours, M) result(R)
    real(dp), intent(in) :: gast_hours, M(3,3)
    real(dp) :: R(3,3)
    R = mat33_mul(rot_z(-gast_hours * TAU / 24.0_dp), M)
  end function

  ! ── Apply IERS polar motion to an ITRS rotation matrix ──
  ! IERS: W = R3(-s') · R2(xp) · R1(yp)  where Ri(θ) = rot_i(-θ)
  ! In our rotation convention: W = rot_y(-xp) · rot_x(-yp)  (dropping s')
  ! To include PM in R_itrs:  R_new = W^(-1) · R_itrs = rot_x(yp) · rot_y(xp) · R_itrs
  ! Matches Skyfield's framelib.py: R = mxm(polar_motion_matrix(), R)
  function apply_polar_motion(R_itrs, xp_as, yp_as) result(R)
    real(dp), intent(in) :: R_itrs(3,3), xp_as, yp_as
    real(dp) :: R(3,3), xp_rad, yp_rad, W_inv(3,3)
    xp_rad = xp_as * ASEC2RAD
    yp_rad = yp_as * ASEC2RAD
    W_inv = mat33_mul(rot_x(yp_rad), rot_y(xp_rad))
    R = mat33_mul(W_inv, R_itrs)
  end function

  ! ── Alt-az rotation ──
  function altaz_rotation(lat_rad, lon_rad, R_itrs) result(R)
    real(dp), intent(in) :: lat_rad, lon_rad, R_itrs(3,3)
    real(dp) :: R(3,3), R_lat(3,3), R_latlon(3,3)
    real(dp) :: temp_row(3)

    R_lat = rot_y(lat_rad)
    ! Reverse row order (Skyfield [::-1] convention)
    temp_row = R_lat(1,:)
    R_lat(1,:) = R_lat(3,:)
    R_lat(3,:) = temp_row

    R_latlon = mat33_mul(R_lat, rot_z(-lon_rad))
    R = mat33_mul(R_latlon, R_itrs)
  end function

  ! ── WGS84 → ITRS (AU) ──
  subroutine wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, pos)
    real(dp), intent(in)  :: lat_deg, lon_deg, elev_m
    real(dp), intent(out) :: pos(3)
    real(dp) :: lat, lon, f, omf2, sinlat, coslat, c, s, xy

    lat = lat_deg * PI / 180.0_dp
    lon = lon_deg * PI / 180.0_dp
    f = 1.0_dp / WGS84_INVF
    omf2 = (1.0_dp - f) ** 2

    sinlat = sin(lat); coslat = cos(lat)
    c = 1.0_dp / sqrt(coslat**2 + sinlat**2 * omf2)
    s = omf2 * c
    xy = (WGS84_RADIUS * c + elev_m) * coslat

    pos(1) = xy * cos(lon)
    pos(2) = xy * sin(lon)
    pos(3) = (WGS84_RADIUS * s + elev_m) * sinlat
    pos = pos / AU_M
  end subroutine

  ! ── ITRS velocity due to Earth rotation ──
  subroutine itrs_velocity_au_per_day(itrs_pos, vel)
    real(dp), intent(in)  :: itrs_pos(3)
    real(dp), intent(out) :: vel(3)
    vel(1) = -itrs_pos(2) * ANGVEL * DAY_S
    vel(2) =  itrs_pos(1) * ANGVEL * DAY_S
    vel(3) = 0.0_dp
  end subroutine

  ! ── Earth position & velocity (SSB, AU) ──
  subroutine earth_position_au(kernel, jd_whole, jd_frac, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    real(dp), intent(out) :: pos(3), vel(3)
    real(dp) :: emb_p(3), emb_v(3), e_p(3), e_v(3)
    call spk_compute_and_diff(kernel, 0, 3, jd_whole, jd_frac, emb_p, emb_v)
    call spk_compute_and_diff(kernel, 3, 399, jd_whole, jd_frac, e_p, e_v)
    pos = (emb_p + e_p) / AU_KM
    vel = (emb_v + e_v) / AU_KM
  end subroutine

  ! ── Sun position & velocity (SSB, AU) ──
  subroutine sun_position_au(kernel, jd_whole, jd_frac, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    real(dp), intent(out) :: pos(3), vel(3)
    call spk_compute_and_diff(kernel, 0, 10, jd_whole, jd_frac, pos, vel)
    pos = pos / AU_KM
    vel = vel / AU_KM
  end subroutine

  ! ── Moon position & velocity (SSB, AU) ──
  subroutine moon_position_au(kernel, jd_whole, jd_frac, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    real(dp), intent(out) :: pos(3), vel(3)
    real(dp) :: emb_p(3), emb_v(3), m_p(3), m_v(3)
    call spk_compute_and_diff(kernel, 0, 3, jd_whole, jd_frac, emb_p, emb_v)
    call spk_compute_and_diff(kernel, 3, 301, jd_whole, jd_frac, m_p, m_v)
    pos = (emb_p + m_p) / AU_KM
    vel = (emb_v + m_v) / AU_KM
  end subroutine

  ! ── Body position at SSB (AU) ──
  subroutine body_ssb_position_au(kernel, jd_whole, jd_frac, target, pos)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    integer, intent(in)   :: target
    real(dp), intent(out) :: pos(3)
    real(dp) :: emb(3), e(3), dummy(3)

    if (target == 10) then
      call spk_compute(kernel, 0, 10, jd_whole, jd_frac, pos)
      pos = pos / AU_KM
    else if (target == 3) then
      call spk_compute(kernel, 0, 3, jd_whole, jd_frac, emb)
      call spk_compute(kernel, 3, 399, jd_whole, jd_frac, e)
      pos = (emb + e) / AU_KM
    else
      call spk_compute(kernel, 0, target, jd_whole, jd_frac, pos)
      pos = pos / AU_KM
    end if
  end subroutine

  ! ── Light-travel-time correction ──
  subroutine correct_light_travel_time(obs_pos, obs_vel, kernel, &
      jd_whole, jd_frac, body_id, astrometric, astro_vel, light_time)
    ! body_id: 1=Sun, 2=Moon
    real(dp), intent(in)  :: obs_pos(3), obs_vel(3)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    integer, intent(in)   :: body_id
    real(dp), intent(out) :: astrometric(3), astro_vel(3), light_time
    real(dp) :: t_pos(3), t_vel(3), dist, lt0, lt, frac2, diff(3)
    integer :: iter

    call get_body_pos_vel(kernel, jd_whole, jd_frac, body_id, t_pos, t_vel)

    diff = t_pos - obs_pos
    dist = vec_len(diff)
    lt0 = 0.0_dp

    do iter = 1, 10
      lt = dist / C_AUDAY
      if (abs(lt - lt0) < 1.0d-12) exit
      frac2 = jd_frac - lt
      call get_body_pos_vel(kernel, jd_whole, frac2, body_id, t_pos, t_vel)
      diff = t_pos - obs_pos
      dist = vec_len(diff)
      lt0 = lt
    end do

    light_time = lt
    astrometric = t_pos - obs_pos
    astro_vel = t_vel - obs_vel
  end subroutine

  ! ── Helper: get body position/velocity by ID ──
  subroutine get_body_pos_vel(kernel, jd_whole, jd_frac, body_id, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, jd_frac
    integer, intent(in)   :: body_id
    real(dp), intent(out) :: pos(3), vel(3)
    if (body_id == 1) then
      call sun_position_au(kernel, jd_whole, jd_frac, pos, vel)
    else
      call moon_position_au(kernel, jd_whole, jd_frac, pos, vel)
    end if
  end subroutine

  ! ── Light-time difference for deflection ──
  function light_time_difference(position, observer_position) result(dlt)
    real(dp), intent(in) :: position(3), observer_position(3)
    real(dp) :: dlt, dis, u1(3)
    dis = vec_len(position)
    if (dis < 1.0d-300) then
      dlt = 0.0_dp
      return
    end if
    u1 = position / dis
    dlt = dot3(u1, observer_position) / C_AUDAY
  end function

  ! ── Single-body deflection ──
  subroutine add_single_deflection(position, observer, deflector, rmass)
    real(dp), intent(inout) :: position(3)
    real(dp), intent(in) :: observer(3), deflector(3), rmass
    real(dp) :: pq(3), pe(3)
    real(dp) :: pmag, qmag, emag
    real(dp) :: phat(3), qhat(3), ehat(3)
    real(dp) :: pdotq, qdote, edotp
    real(dp) :: fac1, fac2
    real(dp), parameter :: AVOID = 1.0d-300

    pq = observer + position - deflector
    pe = observer - deflector

    pmag = vec_len(position)
    qmag = vec_len(pq)
    emag = vec_len(pe)

    phat = position / max(pmag, AVOID)
    qhat = pq / max(qmag, AVOID)
    ehat = pe / max(emag, AVOID)

    pdotq = dot3(phat, qhat)
    qdote = dot3(qhat, ehat)
    edotp = dot3(ehat, phat)

    if (abs(edotp) > 0.99999999999_dp) return

    fac1 = 2.0_dp * GS / (C_SI * C_SI * emag * AU_M * rmass)
    fac2 = 1.0_dp + qdote

    position = position + fac1 * (pdotq * ehat - edotp * qhat) / fac2 * pmag
  end subroutine

  ! ── Limb angle ──
  subroutine compute_limb_angle(position_au, observer_gcrs_au, limb_angle, nadir_angle)
    real(dp), intent(in)  :: position_au(3), observer_gcrs_au(3)
    real(dp), intent(out) :: limb_angle, nadir_angle
    real(dp) :: earth_r_au, disobj, disobs, aprad, zdlim, coszd, zdobj

    earth_r_au = ERAD / AU_M
    disobj = vec_len(position_au)
    disobs = vec_len(observer_gcrs_au)

    aprad = asin(min(earth_r_au / disobs, 1.0_dp))
    zdlim = PI - aprad

    coszd = dot3(position_au, observer_gcrs_au) / (disobj * disobs)
    coszd = max(-1.0_dp, min(1.0_dp, coszd))
    zdobj = acos(coszd)

    limb_angle = (zdlim - zdobj) * (180.0_dp / PI)
    nadir_angle = (PI - zdobj) / aprad
  end subroutine

  ! ── Full deflection (Sun, Jupiter, Saturn + Earth) ──
  subroutine add_deflection(position, obs_bcrs, obs_gcrs, kernel, jd_whole, jd_frac)
    real(dp), intent(inout) :: position(3)
    real(dp), intent(in) :: obs_bcrs(3), obs_gcrs(3)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in) :: jd_whole, jd_frac
    real(dp) :: tlt, bpos(3), gpv(3), dlt, tclose_frac
    real(dp) :: limb_angle, nadir_angle
    real(dp) :: rmasses(3)
    integer :: targets(3), i

    ! Sun, Jupiter, Saturn
    targets = [10, 5, 6]
    rmasses = [RMASS_SUN, RMASS_JUPITER, RMASS_SATURN]

    tlt = vec_len(position) / C_AUDAY

    do i = 1, 3
      call body_ssb_position_au(kernel, jd_whole, jd_frac, targets(i), bpos)
      gpv = bpos - obs_bcrs
      dlt = light_time_difference(position, gpv)

      tclose_frac = jd_frac
      if (dlt > 0.0_dp) tclose_frac = jd_frac - dlt
      if (tlt < dlt) tclose_frac = jd_frac - tlt

      call body_ssb_position_au(kernel, jd_whole, tclose_frac, targets(i), bpos)
      call add_single_deflection(position, obs_bcrs, bpos, rmasses(i))
    end do

    ! Earth deflection
    call compute_limb_angle(position, obs_gcrs, limb_angle, nadir_angle)
    if (nadir_angle >= 0.8_dp) then
      call body_ssb_position_au(kernel, jd_whole, jd_frac, 3, bpos)
      call add_single_deflection(position, obs_bcrs, bpos, RMASS_EARTH)
    end if
  end subroutine

  ! ── Aberration ──
  subroutine add_aberration(position, velocity, light_time)
    real(dp), intent(inout) :: position(3)
    real(dp), intent(in) :: velocity(3), light_time
    real(dp), parameter :: AVOID = 1.0d-300
    real(dp) :: p1mag, vemag, beta, dot_pv, cosd, gammai, p, q, r

    p1mag = light_time * C_AUDAY
    vemag = vec_len(velocity)
    beta = vemag / C_AUDAY
    dot_pv = dot3(position, velocity)

    cosd = dot_pv / (p1mag * vemag + AVOID)
    gammai = sqrt(1.0_dp - beta * beta)
    p = beta * cosd
    q = (1.0_dp + p / (1.0_dp + gammai)) * light_time
    r = 1.0_dp + p

    position = position * gammai
    position = position + q * velocity
    position = position / r
  end subroutine

  ! ── Compute Sun and Moon alt/az/dist at a given UTC time ──
  subroutine compute_altaz_both(kernel, jd_whole, utc_frac, &
      lat_deg, lon_deg, elev_m, delta_t, leap_sec, &
      s_alt, s_az, s_dist, m_alt, m_az, m_dist)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, utc_frac
    real(dp), intent(in)  :: lat_deg, lon_deg, elev_m, delta_t
    integer,  intent(in)  :: leap_sec
    real(dp), intent(out) :: s_alt, s_az, s_dist
    real(dp), intent(out) :: m_alt, m_az, m_dist

    real(dp) :: tt_frac, tdb_frac, ut1_frac
    real(dp) :: jd_tt, jd_tdb
    real(dp) :: Mmat(3,3), d_psi, d_eps, mean_ob
    real(dp) :: gmst_h, gast_h
    real(dp) :: R_it(3,3), RT(3,3), R_aa(3,3)
    real(dp) :: it_pos(3), it_vel(3)
    real(dp) :: o_gcrs(3), o_vgcrs(3)
    real(dp) :: e_pos(3), e_vel(3)
    real(dp) :: o_bcrs(3), o_vbcrs(3)
    real(dp) :: lat_r, lon_r
    real(dp) :: astro(3), astro_v(3), lt
    real(dp) :: aavec(3), dist, alt, az

    ! Time conversions
    call utc_to_tt(jd_whole, utc_frac, leap_sec, tt_frac)
    call tt_to_tdb(jd_whole, tt_frac, tdb_frac)
    call tt_to_ut1(jd_whole, tt_frac, delta_t, ut1_frac)
    jd_tt  = jd_whole + tt_frac
    jd_tdb = jd_whole + tdb_frac

    ! Precession / nutation / sidereal time
    call compute_M(jd_tt, jd_tdb, Mmat, d_psi, d_eps, mean_ob)
    gmst_h = greenwich_mean_sidereal_time(jd_whole, ut1_frac, jd_tdb)
    gast_h = greenwich_apparent_sidereal_time(gmst_h, d_psi, mean_ob, jd_tt)
    R_it = itrs_rotation(gast_h, Mmat)

    ! Observer position (GCRS + BCRS)
    call wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, it_pos)
    call itrs_velocity_au_per_day(it_pos, it_vel)
    RT = mat33_T(R_it)
    o_gcrs  = mat33_vec(RT, it_pos)
    o_vgcrs = mat33_vec(RT, it_vel)
    call earth_position_au(kernel, jd_whole, tdb_frac, e_pos, e_vel)
    o_bcrs  = e_pos + o_gcrs
    o_vbcrs = e_vel + o_vgcrs

    ! Alt-az rotation
    lat_r = lat_deg * PI / 180.0_dp
    lon_r = lon_deg * PI / 180.0_dp
    R_aa = altaz_rotation(lat_r, lon_r, R_it)

    ! ── Sun ──
    call correct_light_travel_time(o_bcrs, o_vbcrs, kernel, &
         jd_whole, tdb_frac, 1, astro, astro_v, lt)
    call add_deflection(astro, o_bcrs, o_gcrs, kernel, jd_whole, tdb_frac)
    call add_aberration(astro, o_vbcrs, lt)
    aavec = mat33_vec(R_aa, astro)
    call to_spherical(aavec, dist, alt, az)
    s_alt  = alt * 180.0_dp / PI
    s_az   = az  * 180.0_dp / PI
    s_dist = dist

    ! ── Moon ──
    call correct_light_travel_time(o_bcrs, o_vbcrs, kernel, &
         jd_whole, tdb_frac, 2, astro, astro_v, lt)
    call add_deflection(astro, o_bcrs, o_gcrs, kernel, jd_whole, tdb_frac)
    call add_aberration(astro, o_vbcrs, lt)
    aavec = mat33_vec(R_aa, astro)
    call to_spherical(aavec, dist, alt, az)
    m_alt  = alt * 180.0_dp / PI
    m_az   = az  * 180.0_dp / PI
    m_dist = dist
  end subroutine

end module astro_mod
