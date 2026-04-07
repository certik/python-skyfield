program artemis_trajectory
  use constants_mod
  use spk_reader_mod
  implicit none

  ! Time conversion constants
  integer, parameter  :: LEAP_SEC = 37           ! UTC leap seconds (2026)
  real(dp), parameter :: TT_OFFSET = (37.0_dp + 32.184_dp) / DAY_S  ! UTC→TT in days

  type(spk_kernel) :: kernel

  integer, parameter  :: MDT_OFFSET = -6 * 3600  ! Mountain Daylight Time = UTC−6

  integer :: n_pts, launch_ts
  integer :: i, u_in, u_out
  integer :: ts_i
  real(dp) :: ox, oy, oz, vel_mag, rr
  real(dp) :: jd_utc, jd_whole, jd_tt_frac, jd_tdb_frac
  real(dp) :: moon_emb(3), earth_emb(3), moon_eq(3), moon_pos(3)
  real(dp) :: earth_dist, moon_dist, met_s
  real(dp) :: dx, dy, dz
  real(dp) :: cos_eps, sin_eps

  ! Max distance tracking
  real(dp) :: max_earth_dist
  integer  :: max_ts, max_met_days, max_met_hrs, max_met_mins
  real(dp) :: max_met_s
  integer  :: mdt_ts, mdt_jd_i, mdt_yr, mdt_mo, mdt_day, mdt_hr, mdt_min, mdt_sec

  ! Precompute obliquity rotation (equatorial → ecliptic)
  cos_eps = cos(OBLIQUITY_RAD)
  sin_eps = sin(OBLIQUITY_RAD)

  ! Read Orion spacecraft data
  print *, "Reading Orion data..."
  open(newunit=u_in, file='artemis_orion.dat', status='old')
  read(u_in, *) n_pts, launch_ts
  print *, "  Records:", n_pts, "  Launch TS:", launch_ts

  ! Open SPK ephemeris
  print *, "Opening de440s.bsp..."
  call spk_open('de440s.bsp', kernel)

  ! Open output file
  open(newunit=u_out, file='artemis_trajectory.dat', status='replace')
  write(u_out, '(A)') '# Artemis II trajectory (Fortran-computed)'
  write(u_out, '(A)') '# Columns: timestamp met_s orion_x orion_y orion_z ' // &
                       'moon_x moon_y moon_z earth_dist_km moon_dist_km velocity_km_s range_rate_km_s'

  max_earth_dist = 0.0_dp
  max_ts = 0

  print *, "Computing trajectory..."
  do i = 1, n_pts
    ! Read one Orion data point
    read(u_in, *) ts_i, ox, oy, oz, vel_mag, rr

    ! Convert Unix timestamp → JD (UTC) → JD (TDB)
    jd_utc = real(ts_i, dp) / DAY_S + UNIX_JD_EPOCH
    jd_whole = floor(jd_utc)
    jd_tt_frac = (jd_utc - jd_whole) + TT_OFFSET
    jd_tdb_frac = jd_tt_frac + tdb_minus_tt(jd_whole, jd_tt_frac) / DAY_S

    ! Compute Moon position relative to Earth (km) in equatorial J2000
    ! Moon_rel_Earth = Moon_rel_EMB − Earth_rel_EMB
    call spk_compute(kernel, 3, 301, jd_whole, jd_tdb_frac, moon_emb)
    call spk_compute(kernel, 3, 399, jd_whole, jd_tdb_frac, earth_emb)
    moon_eq = moon_emb - earth_emb

    ! Rotate equatorial → ecliptic J2000 (Horizons data is ecliptic)
    moon_pos(1) =  moon_eq(1)
    moon_pos(2) =  cos_eps * moon_eq(2) + sin_eps * moon_eq(3)
    moon_pos(3) = -sin_eps * moon_eq(2) + cos_eps * moon_eq(3)

    ! Distances
    earth_dist = sqrt(ox*ox + oy*oy + oz*oz)
    dx = ox - moon_pos(1)
    dy = oy - moon_pos(2)
    dz = oz - moon_pos(3)
    moon_dist = sqrt(dx*dx + dy*dy + dz*dz)

    ! Mission Elapsed Time
    met_s = real(ts_i - launch_ts, dp)

    ! Track maximum Earth distance
    if (earth_dist > max_earth_dist) then
      max_earth_dist = earth_dist
      max_ts = ts_i
    end if

    ! Write output row
    write(u_out, '(I12, 1X, F12.1, 10(1X, ES18.10))') &
        ts_i, met_s, ox, oy, oz, &
        moon_pos(1), moon_pos(2), moon_pos(3), &
        earth_dist, moon_dist, vel_mag, rr
  end do

  close(u_in)
  close(u_out)
  call spk_close(kernel)

  print *, "Wrote", n_pts, " points to artemis_trajectory.dat"

  ! Report maximum Earth distance
  max_met_s = real(max_ts - launch_ts, dp)
  max_met_days = int(max_met_s) / 86400
  max_met_hrs  = mod(int(max_met_s), 86400) / 3600
  max_met_mins = mod(int(max_met_s), 3600) / 60

  ! Convert max_ts (Unix UTC) to Mountain Daylight Time calendar
  mdt_ts = max_ts + MDT_OFFSET
  call unix_to_calendar(mdt_ts, mdt_yr, mdt_mo, mdt_day, mdt_hr, mdt_min, mdt_sec)

  print *, ""
  print *, "═══ Maximum Earth distance ═══"
  write(*, '(A, F12.1, A)') "  Distance:  ", max_earth_dist, " km"
  write(*, '(A, I0, A, I0, A, I2.2, A, I2.2)') &
      "  MET:       T+", max_met_days, "d ", max_met_hrs, "h ", max_met_mins, "m"
  write(*, '(A, I4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A)') &
      "  Time:      ", mdt_yr, "-", mdt_mo, "-", mdt_day, &
      " ", mdt_hr, ":", mdt_min, ":", mdt_sec, " MDT"

contains

  ! Convert Unix timestamp to calendar date/time
  subroutine unix_to_calendar(ts, yr, mo, dy, hr, mn, sc)
    integer, intent(in)  :: ts
    integer, intent(out) :: yr, mo, dy, hr, mn, sc
    integer :: jd, L, N, I_val, J_val, days, secs
    real(dp) :: jd_real

    ! Unix timestamp to Julian Day Number
    days = ts / 86400
    secs = ts - days * 86400
    if (secs < 0) then
      days = days - 1
      secs = secs + 86400
    end if
    jd = days + 2440588  ! JD of Unix day 0 is 2440587.5, noon is +1

    ! Julian Day → calendar (algorithm from Meeus, Astronomical Algorithms)
    L = jd + 68569
    N = 4 * L / 146097
    L = L - (146097 * N + 3) / 4
    I_val = 4000 * (L + 1) / 1461001
    L = L - 1461 * I_val / 4 + 31
    J_val = 80 * L / 2447
    dy = L - 2447 * J_val / 80
    L = J_val / 11
    mo = J_val + 2 - 12 * L
    yr = 100 * (N - 49) + I_val + L

    hr = secs / 3600
    mn = mod(secs, 3600) / 60
    sc = mod(secs, 60)
  end subroutine

  ! TDB − TT correction (USNO Circular 179, eq 2.6)
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

end program artemis_trajectory
