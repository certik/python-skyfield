! ═══════════════════════════════════════════════════════════════════════
!  test_compute_from_scratch.f90
!
!  State-of-the-art apparent-place pipeline using the modern IAU 2006/
!  2000A CIO-based framework with IERS EOP2 corrections:
!
!    [TRS] = RPOM × Rz(ERA) × Q(CIP, s) × [CRS]
!
!  Components:
!    - IAU 2006 precession + IAU 2000A nutation (Fukushima-Williams NPB)
!    - CIP corrected by IERS dX/dY from JPL EOP2
!    - Earth Rotation Angle (ERA, linear in UT1)
!    - Full polar motion matrix RPOM (xp, yp, TIO locator sp)
!    - UT1 from JPL EOP2 (TAI-UT1 + 32.184)
!
!  This is the current IERS standard and yields the highest accuracy
!  against independent implementations (Skyfield/ERFA: < 0.001").
!
!  Note: JPL Horizons uses a DIFFERENT (legacy) pipeline internally —
!  IAU 1976 precession, IAU 1980 nutation (106 terms, GPS-corrected),
!  GMST82 + EQEQ94, and longitude-only polar wobble.  This produces
!  ~0.01–0.5" differences from our results depending on whether their
!  GPS nutation corrections and λ-wobble are replicated.  The modern
!  pipeline here is more accurate than Horizons' legacy framework.
!
!  Checks:
!  1. SKYFIELD/ERFA — computed values must agree with Skyfield (same IAU
!     2006/2000A standard) to < 0.001" for angles and < 1e-13 AU for
!     distances.  Reference values from compute_reference.py.
!  2. HORIZONS — values must agree with JPL Horizons (DE441) to within
!     0.5" for angles (known IAU model difference) and 1e-9 AU for
!     distances.
!  3. ECLIPSE GEOMETRY — Sun-Moon separation confirms eclipse overlap.
!
!  Test case 1: 40°N, 0°E, 0 m — 2025-01-01 12:00 UTC
!  Test case 2: Fredericksburg TX — 2024-04-08 18:35:07 UTC (eclipse max)
!
!  Requires: de441s.bsp, latest_eop2.long, nutation.dat
! ═══════════════════════════════════════════════════════════════════════
program test_compute_from_scratch
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  use eop_mod
  use astro_mod
  use cio_mod
  implicit none

  ! ── Kernel & observer ──────────────────────────────────────────────
  type(spk_kernel) :: kernel
  real(dp) :: lat_deg, lon_deg, elev_m
  integer  :: utc_year, utc_month, utc_day
  integer  :: utc_hour, utc_minute, utc_second
  real(dp) :: delta_t

  ! ── Time ───────────────────────────────────────────────────────────
  real(dp) :: utc_frac, tt_frac, tdb_frac, ut1_frac
  integer  :: jd_int
  real(dp) :: jd_whole, jd_tt, jd_tdb, mjd

  ! ── Orientation ────────────────────────────────────────────────────
  real(dp) :: M(3,3), d_psi, d_eps, mean_ob
  real(dp) :: xp_as, yp_as, tai_ut1_s, dX_mas, dY_mas
  real(dp) :: cip_X, cip_Y, s_cio, sp, era_rad
  real(dp) :: Q(3,3), RPOM(3,3)
  real(dp) :: R_itrs(3,3), RT(3,3), R_altaz(3,3)
  real(dp) :: itrs_pos(3), itrs_vel(3)
  real(dp) :: obs_gcrs(3), obs_vel_gcrs(3)
  real(dp) :: earth_pos(3), earth_vel(3)
  real(dp) :: obs_bcrs_pos(3), obs_bcrs_vel(3)
  real(dp) :: lat_rad, lon_rad

  ! ── Results ────────────────────────────────────────────────────────
  real(dp) :: sun_astro(3), sun_astro_vel(3), sun_lt
  real(dp) :: sun_altaz(3), sun_dist, sun_alt, sun_az
  real(dp) :: sun_ang_diam_as, sun_alt_deg, sun_az_deg

  real(dp) :: moon_astro(3), moon_astro_vel(3), moon_lt
  real(dp) :: moon_altaz(3), moon_dist, moon_alt, moon_az
  real(dp) :: moon_ang_diam_as, moon_alt_deg, moon_az_deg

  real(dp) :: sep_as

  integer :: n_fail

  ! ── Load data files ────────────────────────────────────────────────
  call load_nutation('nutation.dat')
  call load_jpl_eop('latest_eop2.long')
  call spk_open('de441s.bsp', kernel)

  ! ══════════════════════════════════════════════════════════════════
  !  Test 1: 40°N, 0°E — 2025-01-01 12:00 UTC
  ! ══════════════════════════════════════════════════════════════════
  lat_deg = 40.0_dp;  lon_deg = 0.0_dp;  elev_m = 0.0_dp
  utc_year = 2025;  utc_month = 1;  utc_day = 1
  utc_hour = 12;  utc_minute = 0;  utc_second = 0

  jd_int   = julian_day(utc_year, utc_month, utc_day)
  jd_whole = real(jd_int, dp)
  utc_frac = (real(utc_hour,   dp) * 3600.0_dp + &
              real(utc_minute, dp) * 60.0_dp   + &
              real(utc_second, dp)) / DAY_S - 0.5_dp

  ! EOP2: delta_T = 32.184 + TAI-UT1, plus dX/dY and polar motion
  mjd = jd_whole + utc_frac - 2400000.5_dp
  call get_jpl_eop(mjd, xp_as, yp_as, tai_ut1_s, dX_mas, dY_mas)
  delta_t = 32.184_dp + tai_ut1_s

  call utc_to_tt (jd_whole, utc_frac, 37, tt_frac)
  call tt_to_tdb (jd_whole, tt_frac,      tdb_frac)
  call tt_to_ut1 (jd_whole, tt_frac, delta_t, ut1_frac)
  jd_tt  = jd_whole + tt_frac
  jd_tdb = jd_whole + tdb_frac

  ! CIO-based GCRS→ITRS: [TRS] = RPOM × Rz(ERA) × Q × [CRS]
  call compute_npb_fw(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
  cip_X = M(3,1) + (dX_mas / 1000.0_dp) * ASEC2RAD
  cip_Y = M(3,2) + (dY_mas / 1000.0_dp) * ASEC2RAD

  s_cio = compute_cio_s(jd_tt, cip_X, cip_Y)
  Q     = build_cio_matrix(cip_X, cip_Y, s_cio)

  era_rad = earth_rotation_angle(jd_whole, ut1_frac) * TAU

  sp   = tio_locator_sp(jd_tt)
  RPOM = cio_polar_motion(xp_as, yp_as, sp)

  R_itrs = cio_itrs_rotation(Q, era_rad, RPOM)

  call wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, itrs_pos)
  call itrs_velocity_au_per_day(itrs_pos, itrs_vel)
  RT           = mat33_T(R_itrs)
  obs_gcrs     = mat33_vec(RT, itrs_pos)
  obs_vel_gcrs = mat33_vec(RT, itrs_vel)
  call earth_position_au(kernel, jd_whole, tdb_frac, earth_pos, earth_vel)
  obs_bcrs_pos = earth_pos + obs_gcrs
  obs_bcrs_vel = earth_vel + obs_vel_gcrs

  lat_rad = lat_deg * DEG2RAD
  lon_rad = lon_deg * DEG2RAD
  R_altaz = altaz_rotation(lat_rad, lon_rad, R_itrs)

  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 1, sun_astro, sun_astro_vel, sun_lt)
  call add_deflection(sun_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(sun_astro, obs_bcrs_vel, sun_lt)
  sun_altaz       = mat33_vec(R_altaz, sun_astro)
  call to_spherical(sun_altaz, sun_dist, sun_alt, sun_az)
  sun_ang_diam_as = 2.0_dp * asin(SOLAR_RADIUS_KM / (sun_dist * AU_KM)) * RAD2DEG * 3600.0_dp
  sun_alt_deg     = sun_alt * RAD2DEG
  sun_az_deg      = sun_az  * RAD2DEG

  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 2, moon_astro, moon_astro_vel, moon_lt)
  call add_deflection(moon_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(moon_astro, obs_bcrs_vel, moon_lt)
  moon_altaz       = mat33_vec(R_altaz, moon_astro)
  call to_spherical(moon_altaz, moon_dist, moon_alt, moon_az)
  moon_ang_diam_as = 2.0_dp * asin(MOON_RADIUS_KM / (moon_dist * AU_KM)) * RAD2DEG * 3600.0_dp
  moon_alt_deg     = moon_alt * RAD2DEG
  moon_az_deg      = moon_az  * RAD2DEG

  ! ── Assertions: Test 1 ─────────────────────────────────────────────
  n_fail = 0

  print '(A)', '--- Test 1: 40N 0E — 2025-01-01 12:00 UTC ---'

  ! Reference: compute_reference.py (Skyfield + DE441s + EOP2 delta_T & PM)
  ! Both implementations use IAU 2006/2000A; differences arise only from
  ! equinox-vs-CIO framework and minor algorithmic details.
  print '(A)', '=== Skyfield/ERFA reference (IAU 2006/2000A) ==='

  call chk_deg('Sun  alt  vs Skyfield', sun_alt_deg,  27.035994125561917_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_deg('Sun  az   vs Skyfield', sun_az_deg,  179.049495457992748_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_au ('Sun  dist vs Skyfield', sun_dist,      0.983327081438222672_dp, 1.0e-13_dp, n_fail)
  call chk_as ('Sun  diam vs Skyfield', sun_ang_diam_as, 1950.991328191689490_dp, 1.0e-3_dp, n_fail)

  call chk_deg('Moon alt  vs Skyfield', moon_alt_deg, 21.518667068532995_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_deg('Moon az   vs Skyfield', moon_az_deg, 157.820111940143590_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_au ('Moon dist vs Skyfield', moon_dist,     0.00252475127961758732_dp, 1.0e-13_dp, n_fail)
  call chk_as ('Moon diam vs Skyfield', moon_ang_diam_as, 1897.634050540450289_dp, 1.0e-3_dp, n_fail)

  ! Horizons uses a legacy pipeline (IAU 1976/80 + GMST82 + λ-only polar
  ! wobble), so ~0.14" alt / ~0.4" az model-level differences are expected.
  print '(A)', ''
  print '(A)', '=== Horizons accuracy (legacy IAU76/80 pipeline, ~0.5" expected) ==='

  call chk_deg('Sun  alt  vs Horizons', sun_alt_deg,  27.036034_dp,         0.5_dp/3600.0_dp, n_fail)
  call chk_deg('Sun  az   vs Horizons', sun_az_deg,  179.049603_dp,         0.5_dp/3600.0_dp, n_fail)
  call chk_au ('Sun  dist vs Horizons', sun_dist,      0.98332708143732_dp, 1.0e-9_dp,        n_fail)
  call chk_deg('Moon alt  vs Horizons', moon_alt_deg, 21.518703_dp,         0.5_dp/3600.0_dp, n_fail)
  call chk_deg('Moon az   vs Horizons', moon_az_deg, 157.820214_dp,         0.5_dp/3600.0_dp, n_fail)
  call chk_au ('Moon dist vs Horizons', moon_dist,     0.00252475127904_dp, 1.0e-9_dp,        n_fail)
  call chk_as ('Sun  diam vs Horizons', sun_ang_diam_as, 1950.991_dp,       1.0_dp,           n_fail)
  call chk_as ('Moon diam vs Horizons', moon_ang_diam_as, 1897.634_dp,      1.0_dp,           n_fail)

  ! ══════════════════════════════════════════════════════════════════
  !  Test 2: Fredericksburg TX — 2024 April 8 total solar eclipse
  !  Local time 13:35:07 CDT (UTC−5) = 18:35:07 UTC
  ! ══════════════════════════════════════════════════════════════════
  print '(A)', ''
  print '(A)', '--- Test 2: Fredericksburg TX — 2024-04-08 total solar eclipse ---'

  lat_deg = 30.2752011_dp;  lon_deg = -98.8719843_dp;  elev_m = 556.0_dp
  utc_year = 2024;  utc_month = 4;  utc_day = 8
  utc_hour = 18;  utc_minute = 35;  utc_second = 7

  jd_int   = julian_day(utc_year, utc_month, utc_day)
  jd_whole = real(jd_int, dp)
  utc_frac = (real(utc_hour,   dp) * 3600.0_dp + &
              real(utc_minute, dp) * 60.0_dp   + &
              real(utc_second, dp)) / DAY_S - 0.5_dp

  mjd = jd_whole + utc_frac - 2400000.5_dp
  call get_jpl_eop(mjd, xp_as, yp_as, tai_ut1_s, dX_mas, dY_mas)
  delta_t = 32.184_dp + tai_ut1_s

  call utc_to_tt (jd_whole, utc_frac, 37, tt_frac)
  call tt_to_tdb (jd_whole, tt_frac,             tdb_frac)
  call tt_to_ut1 (jd_whole, tt_frac, delta_t,   ut1_frac)
  jd_tt  = jd_whole + tt_frac
  jd_tdb = jd_whole + tdb_frac

  call compute_npb_fw(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
  cip_X = M(3,1) + (dX_mas / 1000.0_dp) * ASEC2RAD
  cip_Y = M(3,2) + (dY_mas / 1000.0_dp) * ASEC2RAD

  s_cio = compute_cio_s(jd_tt, cip_X, cip_Y)
  Q     = build_cio_matrix(cip_X, cip_Y, s_cio)

  era_rad = earth_rotation_angle(jd_whole, ut1_frac) * TAU

  sp   = tio_locator_sp(jd_tt)
  RPOM = cio_polar_motion(xp_as, yp_as, sp)

  R_itrs = cio_itrs_rotation(Q, era_rad, RPOM)

  call wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, itrs_pos)
  call itrs_velocity_au_per_day(itrs_pos, itrs_vel)
  RT           = mat33_T(R_itrs)
  obs_gcrs     = mat33_vec(RT, itrs_pos)
  obs_vel_gcrs = mat33_vec(RT, itrs_vel)
  call earth_position_au(kernel, jd_whole, tdb_frac, earth_pos, earth_vel)
  obs_bcrs_pos = earth_pos + obs_gcrs
  obs_bcrs_vel = earth_vel + obs_vel_gcrs

  lat_rad = lat_deg * DEG2RAD
  lon_rad = lon_deg * DEG2RAD
  R_altaz = altaz_rotation(lat_rad, lon_rad, R_itrs)

  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 1, sun_astro, sun_astro_vel, sun_lt)
  call add_deflection(sun_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(sun_astro, obs_bcrs_vel, sun_lt)
  sun_altaz       = mat33_vec(R_altaz, sun_astro)
  call to_spherical(sun_altaz, sun_dist, sun_alt, sun_az)
  sun_ang_diam_as = 2.0_dp * asin(SOLAR_RADIUS_KM / (sun_dist * AU_KM)) * RAD2DEG * 3600.0_dp
  sun_alt_deg     = sun_alt * RAD2DEG
  sun_az_deg      = sun_az  * RAD2DEG

  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 2, moon_astro, moon_astro_vel, moon_lt)
  call add_deflection(moon_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(moon_astro, obs_bcrs_vel, moon_lt)
  moon_altaz       = mat33_vec(R_altaz, moon_astro)
  call to_spherical(moon_altaz, moon_dist, moon_alt, moon_az)
  moon_ang_diam_as = 2.0_dp * asin(MOON_RADIUS_KM / (moon_dist * AU_KM)) * RAD2DEG * 3600.0_dp
  moon_alt_deg     = moon_alt * RAD2DEG
  moon_az_deg      = moon_az  * RAD2DEG

  ! ── Assertions: Test 2 ─────────────────────────────────────────────

  ! Reference: compute_reference.py (Skyfield + DE441s + EOP2)
  print '(A)', '=== Skyfield/ERFA reference (IAU 2006/2000A) ==='

  call chk_deg('Sun  alt  vs Skyfield', sun_alt_deg,  67.315014549723372_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_deg('Sun  az   vs Skyfield', sun_az_deg,  178.714642811011657_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_au ('Sun  dist vs Skyfield', sun_dist,      1.00147118498815613_dp, 1.0e-13_dp, n_fail)
  call chk_as ('Sun  diam vs Skyfield', sun_ang_diam_as, 1915.644085049606474_dp, 1.0e-3_dp, n_fail)

  call chk_deg('Moon alt  vs Skyfield', moon_alt_deg, 67.316002867115785_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_deg('Moon az   vs Skyfield', moon_az_deg, 178.718239425447138_dp, 0.001_dp/3600.0_dp, n_fail)
  call chk_au ('Moon dist vs Skyfield', moon_dist,     0.00236588022794565010_dp, 1.0e-13_dp, n_fail)
  call chk_as ('Moon diam vs Skyfield', moon_ang_diam_as, 2025.062928456270583_dp, 1.0e-3_dp, n_fail)

  ! Angular separation between Sun and Moon centers (eclipse geometry check)
  ! cos(sep) = sin(alt_s)*sin(alt_m) + cos(alt_s)*cos(alt_m)*cos(az_s - az_m)
  ! At maximum totality the centers should be within a few arcseconds of each other.
  sep_as = acos(min(1.0_dp, &
                    sin(sun_alt) * sin(moon_alt) + &
                    cos(sun_alt) * cos(moon_alt) * cos(sun_az - moon_az))) &
           * RAD2DEG * 3600.0_dp
  print '(A)', ''
  print '(A)', '=== Eclipse geometry check ==='
  print '(A,F8.2,A)', '  Sun-Moon angular separation: ', sep_as, '"'
  print '(A,F8.2,A)', '  Moon angular radius:          ', moon_ang_diam_as / 2.0_dp, '"'
  print '(A,F8.2,A)', '  Sun  angular radius:          ', sun_ang_diam_as  / 2.0_dp, '"'
  call chk_as('Sun-Moon sep < 30"  (eclipse overlap)', sep_as, 0.0_dp, 30.0_dp, n_fail)

  call spk_close(kernel)

  print '(A)', ''
  if (n_fail == 0) then
    print '(A)', 'PASS — all checks passed.'
  else
    print '(I0,A)', n_fail, ' check(s) FAILED.'
    error stop 1
  end if

contains

  subroutine chk_deg(label, got, ref, tol, n_fail)
    character(len=*), intent(in)    :: label
    real(dp),         intent(in)    :: got, ref, tol
    integer,          intent(inout) :: n_fail
    real(dp) :: diff_as
    diff_as = abs(got - ref) * 3600.0_dp   ! arcseconds
    if (abs(got - ref) <= tol) then
      print '(A,A,F10.4,A)', '  PASS  ', label, diff_as, '"'
    else
      print '(A,A,F10.4,A,F10.4,A)', '  FAIL  ', label, diff_as, '"  (tol=', tol*3600.0_dp, '")'
      n_fail = n_fail + 1
    end if
  end subroutine

  subroutine chk_au(label, got, ref, tol, n_fail)
    character(len=*), intent(in)    :: label
    real(dp),         intent(in)    :: got, ref, tol
    integer,          intent(inout) :: n_fail
    real(dp) :: diff
    diff = abs(got - ref)
    if (diff <= tol) then
      print '(A,A,ES10.2,A)', '  PASS  ', label, diff, ' AU'
    else
      print '(A,A,ES10.2,A,ES10.2,A)', '  FAIL  ', label, diff, ' AU  (tol=', tol, ' AU)'
      n_fail = n_fail + 1
    end if
  end subroutine

  subroutine chk_as(label, got, ref, tol, n_fail)
    character(len=*), intent(in)    :: label
    real(dp),         intent(in)    :: got, ref, tol
    integer,          intent(inout) :: n_fail
    real(dp) :: diff
    diff = abs(got - ref)
    if (diff <= tol) then
      print '(A,A,ES10.2,A)', '  PASS  ', label, diff, '"'
    else
      print '(A,A,ES10.2,A,ES10.2,A)', '  FAIL  ', label, diff, '"  (tol=', tol, '")'
      n_fail = n_fail + 1
    end if
  end subroutine

end program test_compute_from_scratch
