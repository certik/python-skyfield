program generate_observations
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  use astro_mod
  implicit none

  type(spk_kernel) :: kernel

  ! ── Observer parameters (edit these) ──
  real(dp) :: lat_deg, lon_deg, elev_m
  real(dp) :: delta_t
  integer  :: leap_sec
  integer  :: year_start, n_days

  ! ── Scan arrays (144 = every 10 min over 24 h) ──
  integer, parameter  :: N_SCAN = 144
  integer, parameter  :: OBS_STEP = 12     ! every 12 scan pts = 2 h
  real(dp) :: sa(0:N_SCAN), saz(0:N_SCAN), sd(0:N_SCAN)  ! Sun
  real(dp) :: ma(0:N_SCAN), maz(0:N_SCAN), md(0:N_SCAN)  ! Moon

  ! ── Noise parameters ──
  real(dp), parameter :: BIAS_HALF  = 0.0_dp  ! no systematic bias (exact obs)
  real(dp), parameter :: SIGMA_BASE = 0.0_dp  ! no random noise (exact obs)
  real(dp), parameter :: SIGMA_HOR  = 0.0_dp  ! no horizon noise (exact obs)
  real(dp) :: bias_sa, bias_saz, bias_ma, bias_maz          ! systematic biases

  ! ── Work variables ──
  integer  :: jd_start, day, i, u_out, n_obs_sun, n_obs_moon
  real(dp) :: jd_day, utc_frac
  real(dp) :: s_alt, s_az, s_dist, m_alt, m_az, m_dist
  real(dp) :: cross_frac, c_alt, c_az, c_dist
  real(dp) :: alt_noisy, az_noisy
  real(dp) :: r_val
  integer  :: n_seed
  integer, allocatable :: seed(:)

  character(len=256) :: bsp_file

  ! ── Observer location (Ondrej's backyard — edit as needed) ──
  lat_deg    = 40.0_dp
  lon_deg    = 0.0_dp
  elev_m     = 0.0_dp
  delta_t    = 69.15_dp
  leap_sec   = 37
  year_start = 2025
  n_days     = 365

  ! ── Load ephemeris data ──
  call get_command_argument(1, bsp_file)
  if (len_trim(bsp_file) == 0) bsp_file = 'de440s.bsp'
  call load_nutation('nutation.dat')
  call spk_open(trim(bsp_file), kernel)

  ! ── Initialise reproducible RNG ──
  call random_seed(size=n_seed)
  allocate(seed(n_seed))
  do i = 1, n_seed
    seed(i) = 42 + i * 137
  end do
  call random_seed(put=seed)
  deallocate(seed)

  ! ── Draw systematic biases ──
  call random_number(r_val); bias_sa  = (r_val - 0.5_dp) * 2.0_dp * BIAS_HALF
  call random_number(r_val); bias_saz = (r_val - 0.5_dp) * 2.0_dp * BIAS_HALF
  call random_number(r_val); bias_ma  = (r_val - 0.5_dp) * 2.0_dp * BIAS_HALF
  call random_number(r_val); bias_maz = (r_val - 0.5_dp) * 2.0_dp * BIAS_HALF

  ! ── Open output file ──
  open(newunit=u_out, file='observations.dat', &
       status='replace', action='write')
  write(u_out, '(A)')       '# Synthetic backyard observations for Kalman-filter input'
  write(u_out, '(A,F9.4,A,F9.4,A,F6.1)') &
       '# Observer  lat=', lat_deg, '  lon=', lon_deg, '  elev=', elev_m
  write(u_out, '(A,I4,A,I4)') '# Span ', year_start, '-01-01  days=', n_days
  write(u_out, '(A)')       '# Noise: systematic bias + Gaussian random'
  write(u_out, '(A)')       '# Biases hidden (realistic unknown calibration error)'
  write(u_out, '(A)')       '# Average error < 1 arcmin from true position'
  write(u_out, '(A)')       '#'
  write(u_out, '(A)')       '# Columns:'
  write(u_out, '(A)')       '#  1  JD_UTC          Julian date (UTC)'
  write(u_out, '(A)')       '#  2  body             S=Sun  M=Moon'
  write(u_out, '(A)')       '#  3  event            R=rise  S=set  I=intermediate'
  write(u_out, '(A)')       '#  4  alt_obs (deg)    measured altitude'
  write(u_out, '(A)')       '#  5  az_obs  (deg)    measured azimuth (N=0 E=90)'
  write(u_out, '(A)')       '#  6  alt_true (deg)   true altitude'
  write(u_out, '(A)')       '#  7  az_true  (deg)   true azimuth'
  write(u_out, '(A)')       '#  8  dist     (AU)    geocentric distance'

  n_obs_sun  = 0
  n_obs_moon = 0

  ! ── Header to stdout ──
  print '(A)',        '========================================================'
  print '(A)',        'Generating synthetic observations'
  print '(A,F8.4,A,F8.4)', '  lat=', lat_deg, '  lon=', lon_deg
  print '(A,I4,A,I4)',     '  year=', year_start, '  days=', n_days
  print '(A)',        '========================================================'

  jd_start = julian_day(year_start, 1, 1)

  ! ════════════════════════════════════════════════════════════════════
  !  Day loop
  ! ════════════════════════════════════════════════════════════════════
  do day = 0, n_days - 1
    jd_day = real(jd_start + day, dp)

    ! ── Coarse scan: every 10 min ──
    do i = 0, N_SCAN
      utc_frac = real(i, dp) / real(N_SCAN, dp) - 0.5_dp
      call compute_altaz_both(kernel, jd_day, utc_frac, &
           lat_deg, lon_deg, elev_m, delta_t, leap_sec, &
           s_alt, s_az, s_dist, m_alt, m_az, m_dist)
      sa(i)  = s_alt;  saz(i) = s_az;  sd(i) = s_dist
      ma(i)  = m_alt;  maz(i) = m_az;  md(i) = m_dist
    end do

    ! ── Sun: find rise/set crossings ──
    do i = 1, N_SCAN
      ! Sunrise
      if (sa(i-1) < 0.0_dp .and. sa(i) >= 0.0_dp) then
        call bisect_body(kernel, jd_day, &
             real(i-1, dp)/real(N_SCAN, dp) - 0.5_dp, &
             real(i,   dp)/real(N_SCAN, dp) - 0.5_dp, &
             lat_deg, lon_deg, elev_m, delta_t, leap_sec, 1, &
             cross_frac, c_alt, c_az, c_dist)
        call apply_noise(c_alt, c_az, bias_sa, bias_saz, alt_noisy, az_noisy)
        call write_obs(u_out, jd_day + cross_frac, 'S', 'R', &
             alt_noisy, az_noisy, c_alt, c_az, c_dist)
        n_obs_sun = n_obs_sun + 1
      end if
      ! Sunset
      if (sa(i-1) >= 0.0_dp .and. sa(i) < 0.0_dp) then
        call bisect_body(kernel, jd_day, &
             real(i-1, dp)/real(N_SCAN, dp) - 0.5_dp, &
             real(i,   dp)/real(N_SCAN, dp) - 0.5_dp, &
             lat_deg, lon_deg, elev_m, delta_t, leap_sec, 1, &
             cross_frac, c_alt, c_az, c_dist)
        call apply_noise(c_alt, c_az, bias_sa, bias_saz, alt_noisy, az_noisy)
        call write_obs(u_out, jd_day + cross_frac, 'S', 'S', &
             alt_noisy, az_noisy, c_alt, c_az, c_dist)
        n_obs_sun = n_obs_sun + 1
      end if
    end do

    ! Sun intermediates (every 2 h while above horizon)
    do i = 0, N_SCAN, OBS_STEP
      if (sa(i) > 0.0_dp) then
        call apply_noise(sa(i), saz(i), bias_sa, bias_saz, alt_noisy, az_noisy)
        utc_frac = real(i, dp) / real(N_SCAN, dp) - 0.5_dp
        call write_obs(u_out, jd_day + utc_frac, 'S', 'I', &
             alt_noisy, az_noisy, sa(i), saz(i), sd(i))
        n_obs_sun = n_obs_sun + 1
      end if
    end do

    ! ── Moon: find rise/set crossings ──
    do i = 1, N_SCAN
      ! Moonrise
      if (ma(i-1) < 0.0_dp .and. ma(i) >= 0.0_dp) then
        call bisect_body(kernel, jd_day, &
             real(i-1, dp)/real(N_SCAN, dp) - 0.5_dp, &
             real(i,   dp)/real(N_SCAN, dp) - 0.5_dp, &
             lat_deg, lon_deg, elev_m, delta_t, leap_sec, 2, &
             cross_frac, c_alt, c_az, c_dist)
        call apply_noise(c_alt, c_az, bias_ma, bias_maz, alt_noisy, az_noisy)
        call write_obs(u_out, jd_day + cross_frac, 'M', 'R', &
             alt_noisy, az_noisy, c_alt, c_az, c_dist)
        n_obs_moon = n_obs_moon + 1
      end if
      ! Moonset
      if (ma(i-1) >= 0.0_dp .and. ma(i) < 0.0_dp) then
        call bisect_body(kernel, jd_day, &
             real(i-1, dp)/real(N_SCAN, dp) - 0.5_dp, &
             real(i,   dp)/real(N_SCAN, dp) - 0.5_dp, &
             lat_deg, lon_deg, elev_m, delta_t, leap_sec, 2, &
             cross_frac, c_alt, c_az, c_dist)
        call apply_noise(c_alt, c_az, bias_ma, bias_maz, alt_noisy, az_noisy)
        call write_obs(u_out, jd_day + cross_frac, 'M', 'S', &
             alt_noisy, az_noisy, c_alt, c_az, c_dist)
        n_obs_moon = n_obs_moon + 1
      end if
    end do

    ! Moon intermediates (every 2 h while above horizon)
    do i = 0, N_SCAN, OBS_STEP
      if (ma(i) > 0.0_dp) then
        call apply_noise(ma(i), maz(i), bias_ma, bias_maz, alt_noisy, az_noisy)
        utc_frac = real(i, dp) / real(N_SCAN, dp) - 0.5_dp
        call write_obs(u_out, jd_day + utc_frac, 'M', 'I', &
             alt_noisy, az_noisy, ma(i), maz(i), md(i))
        n_obs_moon = n_obs_moon + 1
      end if
    end do

    ! Progress every 30 days
    if (mod(day, 30) == 0) then
      print '(A,I4,A,I6,A,I6)', '  day ', day, &
            '  sun_obs=', n_obs_sun, '  moon_obs=', n_obs_moon
    end if
  end do

  close(u_out)
  call spk_close(kernel)

  ! ── Summary ──
  print '(/,A)',      '========================================================'
  print '(A)',        'Done.'
  print '(A,I7)',     '  Sun observations:  ', n_obs_sun
  print '(A,I7)',     '  Moon observations: ', n_obs_moon
  print '(A,I7)',     '  Total:             ', n_obs_sun + n_obs_moon
  print '(A)',        '  Output: observations.dat'
  print '(A)',        '========================================================'
  print '(A)',        'Systematic biases (arcsec) — hidden from observer:'
  print '(A,F8.2)',   '  Sun  alt bias: ', bias_sa  * 3600.0_dp
  print '(A,F8.2)',   '  Sun  az  bias: ', bias_saz * 3600.0_dp
  print '(A,F8.2)',   '  Moon alt bias: ', bias_ma  * 3600.0_dp
  print '(A,F8.2)',   '  Moon az  bias: ', bias_maz * 3600.0_dp

contains

  ! ── Bisection to find horizon crossing for one body ──
  subroutine bisect_body(kernel, jd_whole, frac_lo, frac_hi, &
      lat_d, lon_d, elev, dt, lsec, body_idx, &
      frac_out, alt_out, az_out, dist_out)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in)  :: jd_whole, frac_lo, frac_hi
    real(dp), intent(in)  :: lat_d, lon_d, elev, dt
    integer,  intent(in)  :: lsec, body_idx
    real(dp), intent(out) :: frac_out, alt_out, az_out, dist_out

    real(dp) :: lo, hi, mid
    real(dp) :: s_a, s_z, s_d, m_a, m_z, m_d
    real(dp) :: alt_lo, alt_mid
    integer  :: iter

    lo = frac_lo;  hi = frac_hi

    ! altitude at lo
    call compute_altaz_both(kernel, jd_whole, lo, &
         lat_d, lon_d, elev, dt, lsec, s_a, s_z, s_d, m_a, m_z, m_d)
    if (body_idx == 1) then; alt_lo = s_a; else; alt_lo = m_a; end if

    do iter = 1, 20
      mid = (lo + hi) * 0.5_dp
      call compute_altaz_both(kernel, jd_whole, mid, &
           lat_d, lon_d, elev, dt, lsec, s_a, s_z, s_d, m_a, m_z, m_d)
      if (body_idx == 1) then; alt_mid = s_a; else; alt_mid = m_a; end if

      if (alt_lo * alt_mid < 0.0_dp) then
        hi = mid
      else
        lo = mid
        alt_lo = alt_mid
      end if
    end do

    frac_out = (lo + hi) * 0.5_dp
    ! Final evaluation at crossing
    call compute_altaz_both(kernel, jd_whole, frac_out, &
         lat_d, lon_d, elev, dt, lsec, s_a, s_z, s_d, m_a, m_z, m_d)
    if (body_idx == 1) then
      alt_out = s_a; az_out = s_z; dist_out = s_d
    else
      alt_out = m_a; az_out = m_z; dist_out = m_d
    end if
  end subroutine

  ! ── Add systematic bias + random noise ──
  subroutine apply_noise(alt_true, az_true, b_alt, b_az, alt_out, az_out)
    real(dp), intent(in)  :: alt_true, az_true, b_alt, b_az
    real(dp), intent(out) :: alt_out, az_out
    real(dp) :: sigma, u1, u2, g1, g2

    ! Larger noise near the horizon
    if (abs(alt_true) < 5.0_dp) then
      sigma = sqrt(SIGMA_BASE**2 + SIGMA_HOR**2)
    else
      sigma = SIGMA_BASE
    end if

    ! Box-Muller for two independent Gaussian samples
    call random_number(u1); call random_number(u2)
    u1 = max(u1, 1.0d-300)
    g1 = sqrt(-2.0_dp * log(u1)) * cos(TAU * u2)
    g2 = sqrt(-2.0_dp * log(u1)) * sin(TAU * u2)

    alt_out = alt_true + b_alt + sigma * g1
    az_out  = az_true  + b_az  + sigma * g2
    ! Keep azimuth in [0, 360)
    az_out = mod(az_out, 360.0_dp)
    if (az_out < 0.0_dp) az_out = az_out + 360.0_dp
  end subroutine

  ! ── Write one observation line ──
  subroutine write_obs(u, jd, body, event, alt_m, az_m, alt_t, az_t, dist)
    integer,      intent(in) :: u
    real(dp),     intent(in) :: jd, alt_m, az_m, alt_t, az_t, dist
    character(1), intent(in) :: body, event
    write(u, '(F15.6,2X,A1,2X,A1,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,ES14.7)') &
         jd, body, event, alt_m, az_m, alt_t, az_t, dist
  end subroutine

end program generate_observations
