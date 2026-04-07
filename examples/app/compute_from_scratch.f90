program compute_from_scratch
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  use astro_mod
  implicit none

  type(spk_kernel) :: kernel
  real(dp) :: lat_deg, lon_deg, elev_m
  integer  :: utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_second
  real(dp) :: delta_t
  integer  :: leap_sec
  real(dp) :: utc_frac, tt_frac, tdb_frac, ut1_frac
  integer  :: jd_int
  real(dp) :: jd_whole, jd_tt, jd_tdb
  real(dp) :: M(3,3), d_psi, d_eps, mean_ob
  real(dp) :: gmst_h, gast_h
  real(dp) :: R_itrs(3,3), RT(3,3), R_altaz(3,3)
  real(dp) :: itrs_pos(3), itrs_vel(3)
  real(dp) :: obs_gcrs(3), obs_vel_gcrs(3)
  real(dp) :: earth_pos(3), earth_vel(3)
  real(dp) :: obs_bcrs_pos(3), obs_bcrs_vel(3)
  real(dp) :: lat_rad, lon_rad

  ! Sun
  real(dp) :: sun_astro(3), sun_astro_vel(3), sun_lt
  real(dp) :: sun_altaz(3), sun_dist, sun_alt, sun_az
  real(dp) :: sun_ang_diam, sun_alt_deg, sun_az_deg, sun_ang_diam_as

  ! Moon
  real(dp) :: moon_astro(3), moon_astro_vel(3), moon_lt
  real(dp) :: moon_altaz(3), moon_dist, moon_alt, moon_az
  real(dp) :: moon_ang_diam, moon_alt_deg, moon_az_deg, moon_ang_diam_as

  ! Horizons reference
  real(dp) :: diff_val
  character(len=256) :: bsp_file

  ! Load data files
  call get_command_argument(1, bsp_file)
  if (len_trim(bsp_file) == 0) bsp_file = 'de440s.bsp'
  call load_nutation('nutation.dat')
  call spk_open(trim(bsp_file), kernel)

  ! ── Test case parameters ──
  lat_deg = 40.0_dp
  lon_deg = 0.0_dp
  elev_m  = 0.0_dp
  utc_year = 2025; utc_month = 1; utc_day = 1
  utc_hour = 12; utc_minute = 0; utc_second = 0
  delta_t = 69.14980035_dp
  leap_sec = 37

  print '(A)', '================================================================='
  print '(A)', '40 N, Greenwich -- 2025 January 1, 12:00 UTC'
  print '(A)', 'SPK backend: Fortran (pure)'
  print '(A,A)', 'SPK file: ', trim(bsp_file)
  print '(A)', '================================================================='

  ! ── Time conversions ──
  jd_int = julian_day(utc_year, utc_month, utc_day)
  jd_whole = real(jd_int, dp)
  utc_frac = (real(utc_hour, dp) * 3600.0_dp + &
              real(utc_minute, dp) * 60.0_dp + &
              real(utc_second, dp)) / DAY_S - 0.5_dp

  call utc_to_tt(jd_whole, utc_frac, leap_sec, tt_frac)
  call tt_to_tdb(jd_whole, tt_frac, tdb_frac)
  call tt_to_ut1(jd_whole, tt_frac, delta_t, ut1_frac)

  jd_tt  = jd_whole + tt_frac
  jd_tdb = jd_whole + tdb_frac

  print '(A,F24.15)', '  JD (TT):  ', jd_tt
  print '(A,F24.15)', '  JD (TDB): ', jd_tdb
  print '(A,F12.5,A)', '  delta_T:  ', delta_t, ' s'

  ! ── Precession, nutation, sidereal time ──
  call compute_M(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
  gmst_h = greenwich_mean_sidereal_time(jd_whole, ut1_frac, jd_tdb)
  gast_h = greenwich_apparent_sidereal_time(gmst_h, d_psi, mean_ob, jd_tt)
  R_itrs = itrs_rotation(gast_h, M)

  print '(A,F20.15,A)', '  GMST: ', gmst_h, ' h'
  print '(A,F20.15,A)', '  GAST: ', gast_h, ' h'

  ! ── Observer position ──
  call wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, itrs_pos)
  call itrs_velocity_au_per_day(itrs_pos, itrs_vel)

  RT = mat33_T(R_itrs)
  obs_gcrs = mat33_vec(RT, itrs_pos)
  obs_vel_gcrs = mat33_vec(RT, itrs_vel)

  ! ── Earth barycentric position ──
  call earth_position_au(kernel, jd_whole, tdb_frac, earth_pos, earth_vel)

  obs_bcrs_pos = earth_pos + obs_gcrs
  obs_bcrs_vel = earth_vel + obs_vel_gcrs

  print '(A,3ES20.9)', '  Observer BCRS:', obs_bcrs_pos
  print '(A,3ES20.9)', '  Observer GCRS:', obs_gcrs

  ! ── Alt-az rotation ──
  lat_rad = lat_deg * PI / 180.0_dp
  lon_rad = lon_deg * PI / 180.0_dp
  R_altaz = altaz_rotation(lat_rad, lon_rad, R_itrs)

  ! ═══════════════════════════════════════════════════════════════════
  !  Sun
  ! ═══════════════════════════════════════════════════════════════════
  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 1, sun_astro, sun_astro_vel, sun_lt)

  print '(/,A,3ES20.9)', '  Sun astrometric:', sun_astro
  print '(A,F20.15,A)', '  Sun light time:  ', sun_lt, ' days'

  ! Deflection + aberration
  call add_deflection(sun_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(sun_astro, obs_bcrs_vel, sun_lt)

  ! Alt-az
  sun_altaz = mat33_vec(R_altaz, sun_astro)
  call to_spherical(sun_altaz, sun_dist, sun_alt, sun_az)
  sun_ang_diam = 2.0_dp * asin(SOLAR_RADIUS_KM / (sun_dist * AU_KM))
  sun_alt_deg = sun_alt * 180.0_dp / PI
  sun_az_deg  = sun_az * 180.0_dp / PI
  sun_ang_diam_as = sun_ang_diam * 180.0_dp / PI * 3600.0_dp

  ! ═══════════════════════════════════════════════════════════════════
  !  Moon
  ! ═══════════════════════════════════════════════════════════════════
  call correct_light_travel_time(obs_bcrs_pos, obs_bcrs_vel, kernel, &
       jd_whole, tdb_frac, 2, moon_astro, moon_astro_vel, moon_lt)

  print '(A,3ES20.9)', '  Moon astrometric:', moon_astro
  print '(A,ES25.15,A)', '  Moon light time:  ', moon_lt, ' days'

  ! Deflection + aberration
  call add_deflection(moon_astro, obs_bcrs_pos, obs_gcrs, kernel, jd_whole, tdb_frac)
  call add_aberration(moon_astro, obs_bcrs_vel, moon_lt)

  ! Alt-az
  moon_altaz = mat33_vec(R_altaz, moon_astro)
  call to_spherical(moon_altaz, moon_dist, moon_alt, moon_az)
  moon_ang_diam = 2.0_dp * asin(MOON_RADIUS_KM / (moon_dist * AU_KM))
  moon_alt_deg = moon_alt * 180.0_dp / PI
  moon_az_deg  = moon_az * 180.0_dp / PI
  moon_ang_diam_as = moon_ang_diam * 180.0_dp / PI * 3600.0_dp

  ! ═══════════════════════════════════════════════════════════════════
  !  Print results
  ! ═══════════════════════════════════════════════════════════════════
  print '(/,A)', 'Sun:'
  print '(A,F13.6,A)', '  Altitude:  ', sun_alt_deg, ' deg'
  print '(A,F13.6,A)', '  Azimuth:   ', sun_az_deg, ' deg'
  print '(A,F10.3,A)',  '  Ang-diam:  ', sun_ang_diam_as, '"'
  print '(A,F19.14,A)', '  delta:     ', sun_dist, ' AU'

  print '(A)', 'Moon:'
  print '(A,F13.6,A)', '  Altitude:  ', moon_alt_deg, ' deg'
  print '(A,F13.6,A)', '  Azimuth:   ', moon_az_deg, ' deg'
  print '(A,F10.3,A)',  '  Ang-diam:  ', moon_ang_diam_as, '"'
  print '(A,F19.14,A)', '  delta:     ', moon_dist, ' AU'

  ! ═══════════════════════════════════════════════════════════════════
  !  Comparison with JPL Horizons
  ! ═══════════════════════════════════════════════════════════════════
  print '(/,A)', '================================================================='
  print '(A)',   'Comparison with JPL Horizons (DE441)'
  print '(A)',   '================================================================='

  diff_val = sun_alt_deg - 27.036034_dp
  print '(A,F13.6,A,F13.6,A,F8.3,A)', '  Sun alt:      ', sun_alt_deg, &
    '  Horizons: ', 27.036034_dp, '  diff: ', abs(diff_val)*3600.0_dp, '"'
  diff_val = sun_az_deg - 179.049603_dp
  print '(A,F13.6,A,F13.6,A,F8.3,A)', '  Sun az:       ', sun_az_deg, &
    '  Horizons: ', 179.049603_dp, '  diff: ', abs(diff_val)*3600.0_dp, '"'
  diff_val = sun_ang_diam_as - 1950.991_dp
  print '(A,F13.3,A,F13.3,A,F8.3,A)', '  Sun Ang-diam: ', sun_ang_diam_as, &
    '  Horizons: ', 1950.991_dp, '  diff: ', diff_val, '"'
  diff_val = sun_dist - 0.98332708143732_dp
  print '(A,F19.14,A,F19.14,A,ES11.2,A)', '  Sun delta:    ', sun_dist, &
    '  Horizons: ', 0.98332708143732_dp, '  diff: ', diff_val, ' AU'

  diff_val = moon_alt_deg - 21.518703_dp
  print '(A,F13.6,A,F13.6,A,F8.3,A)', '  Moon alt:     ', moon_alt_deg, &
    '  Horizons: ', 21.518703_dp, '  diff: ', abs(diff_val)*3600.0_dp, '"'
  diff_val = moon_az_deg - 157.820214_dp
  print '(A,F13.6,A,F13.6,A,F8.3,A)', '  Moon az:      ', moon_az_deg, &
    '  Horizons: ', 157.820214_dp, '  diff: ', abs(diff_val)*3600.0_dp, '"'
  diff_val = moon_ang_diam_as - 1897.634_dp
  print '(A,F13.3,A,F13.3,A,F8.3,A)', '  Moon Ang-diam:', moon_ang_diam_as, &
    '  Horizons: ', 1897.634_dp, '  diff: ', diff_val, '"'
  diff_val = moon_dist - 0.00252475127904_dp
  print '(A,F19.14,A,F19.14,A,ES11.2,A)', '  Moon delta:   ', moon_dist, &
    '  Horizons: ', 0.00252475127904_dp, '  diff: ', diff_val, ' AU'

  print '(/,A)', '  Note: Horizons uses DE441 ephemeris + its own EOP data.'
  print '(A)',   '  We use DE440s + a single delta_T value. Small alt/az'
  print '(A)',   '  differences (~0.1") are expected from these sources.'

  call spk_close(kernel)

end program compute_from_scratch
