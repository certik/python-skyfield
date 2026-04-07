! ═══════════════════════════════════════════════════════════════════════
!  test_de441_horizons.f90
!
!  Verify that our pipeline matches JPL Horizons when using the SAME
!  ephemeris (DE441).  Horizons uses DE441 internally, so with
!  de441_part-2.bsp the only remaining differences are:
!    - Earth orientation (we use a fixed delta_T; Horizons uses EOP)
!    - Nutation model details
!
!  Distances and angular diameters should agree to machine precision.
!  Angles (alt/az) may still differ by ~0.4" due to Earth orientation.
!
!  Requires: de441_part-2.bsp  (run download_data.sh to obtain it)
!  If the file is absent the test prints SKIP and exits successfully.
! ═══════════════════════════════════════════════════════════════════════
program test_de441_horizons
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  use astro_mod
  implicit none

  type(spk_kernel) :: kernel
  real(dp) :: lat_deg, lon_deg, elev_m
  integer  :: utc_year, utc_month, utc_day
  integer  :: utc_hour, utc_minute, utc_second
  real(dp) :: delta_t
  integer  :: leap_sec

  real(dp) :: utc_frac, tt_frac, tdb_frac, ut1_frac
  integer  :: jd_int
  real(dp) :: jd_whole, jd_tt, jd_tdb

  real(dp) :: M(3,3), d_psi, d_eps, mean_ob, gmst_h, gast_h
  real(dp) :: R_itrs(3,3), RT(3,3), R_altaz(3,3)
  real(dp) :: itrs_pos(3), itrs_vel(3)
  real(dp) :: obs_gcrs(3), obs_vel_gcrs(3)
  real(dp) :: earth_pos(3), earth_vel(3)
  real(dp) :: obs_bcrs_pos(3), obs_bcrs_vel(3)
  real(dp) :: lat_rad, lon_rad

  real(dp) :: sun_astro(3), sun_astro_vel(3), sun_lt
  real(dp) :: sun_altaz(3), sun_dist, sun_alt, sun_az
  real(dp) :: sun_ang_diam_as, sun_alt_deg, sun_az_deg

  real(dp) :: moon_astro(3), moon_astro_vel(3), moon_lt
  real(dp) :: moon_altaz(3), moon_dist, moon_alt, moon_az
  real(dp) :: moon_ang_diam_as, moon_alt_deg, moon_az_deg

  integer :: n_fail
  logical :: file_exists

  ! ── Check for DE441s ─────────────────────────────────────────────────
  inquire(file='de441s.bsp', exist=file_exists)
  if (.not. file_exists) then
    ! Fall back to the full part-2 file
    inquire(file='de441_part-2.bsp', exist=file_exists)
    if (.not. file_exists) then
      print '(A)', 'SKIP — de441s.bsp not found (run: fpm run create_de441s)'
      stop
    end if
    call load_nutation('nutation.dat')
    call spk_open('de441_part-2.bsp', kernel)
  else
    call load_nutation('nutation.dat')
    call spk_open('de441s.bsp', kernel)
  end if

  ! ══════════════════════════════════════════════════════════════════
  !  40°N, 0°E, 0 m — 2025-01-01 12:00 UTC
  !
  !  Horizons reference (DE441, airless, same observer):
  !    Sun  alt  27.036034°   az  179.049603°   dist 0.98332708143732 AU
  !         ang-diam 1950.991"
  !    Moon alt  21.518703°   az  157.820214°   dist 0.00252475127904 AU
  !         ang-diam 1897.634"
  ! ══════════════════════════════════════════════════════════════════
  n_fail = 0

  print '(A)', '--- DE441 Horizons check: 40N 0E — 2025-01-01 12:00 UTC ---'

  lat_deg = 40.0_dp;  lon_deg = 0.0_dp;  elev_m = 0.0_dp
  utc_year = 2025;  utc_month = 1;  utc_day = 1
  utc_hour = 12;  utc_minute = 0;  utc_second = 0
  delta_t  = 69.14980035_dp
  leap_sec = 37

  jd_int   = julian_day(utc_year, utc_month, utc_day)
  jd_whole = real(jd_int, dp)
  utc_frac = (real(utc_hour,   dp) * 3600.0_dp + &
              real(utc_minute, dp) * 60.0_dp   + &
              real(utc_second, dp)) / DAY_S - 0.5_dp

  call utc_to_tt (jd_whole, utc_frac, leap_sec, tt_frac)
  call tt_to_tdb (jd_whole, tt_frac,             tdb_frac)
  call tt_to_ut1 (jd_whole, tt_frac, delta_t,   ut1_frac)
  jd_tt  = jd_whole + tt_frac
  jd_tdb = jd_whole + tdb_frac

  call compute_M(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
  gmst_h = greenwich_mean_sidereal_time(jd_whole, ut1_frac, jd_tdb)
  gast_h = greenwich_apparent_sidereal_time(gmst_h, d_psi, mean_ob, jd_tt)
  R_itrs = itrs_rotation(gast_h, M)

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

  ! ── Assertions ──────────────────────────────────────────────────────
  ! Distances depend only on the ephemeris → should match Horizons
  ! to ~1e-11 AU or better (light-travel-time iteration limit).
  ! Angular diameters derive from distance → same tight tolerance.
  ! Alt/az still differ by ~0.4" due to Earth orientation (delta_T vs EOP).
  print '(A)', '=== Horizons vs DE441 ==='

  call chk_deg('Sun  alt  vs Horizons', sun_alt_deg,  27.036034_dp,         1.0_dp/3600.0_dp, n_fail)
  call chk_deg('Sun  az   vs Horizons', sun_az_deg,  179.049603_dp,         1.0_dp/3600.0_dp, n_fail)
  call chk_au ('Sun  dist vs Horizons', sun_dist,      0.98332708143732_dp, 5.0e-11_dp,       n_fail)
  call chk_as ('Sun  diam vs Horizons', sun_ang_diam_as, 1950.991_dp,       1.0e-2_dp,        n_fail)

  call chk_deg('Moon alt  vs Horizons', moon_alt_deg, 21.518703_dp,         1.0_dp/3600.0_dp, n_fail)
  call chk_deg('Moon az   vs Horizons', moon_az_deg, 157.820214_dp,         1.0_dp/3600.0_dp, n_fail)
  call chk_au ('Moon dist vs Horizons', moon_dist,     0.00252475127904_dp, 1.0e-11_dp,       n_fail)
  call chk_as ('Moon diam vs Horizons', moon_ang_diam_as, 1897.634_dp,      1.0e-2_dp,        n_fail)

  call spk_close(kernel)

  print '(A)', ''
  if (n_fail == 0) then
    print '(A)', 'PASS — all DE441 Horizons checks passed.'
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
    diff_as = abs(got - ref) * 3600.0_dp
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

end program test_de441_horizons
