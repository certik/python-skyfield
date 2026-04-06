! ═══════════════════════════════════════════════════════════════════════
!  artemis_trajectory.f90
!
!  Computes Artemis II mission trajectory data for plotting.
!  Reads Orion spacecraft positions (from JPL Horizons via JSON→text)
!  and computes Moon position from de440s.bsp ephemeris.
!
!  Input:  artemis_orion.dat  (Orion positions, from convert_artemis_data.py)
!          de440s.bsp         (JPL planetary ephemeris)
!  Output: artemis_trajectory.dat  (combined trajectory data for plotting)
!
!  Compile:  lfortran artemis_trajectory.f90 -o artemis_trajectory
!  Run:      ./artemis_trajectory
! ═══════════════════════════════════════════════════════════════════════


! ─────────────────────────────────────────────────────────────────────
!  Constants
! ─────────────────────────────────────────────────────────────────────
module constants_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI      = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: T0      = 2451545.0_dp           ! J2000.0 epoch (JD)
  real(dp), parameter :: DAY_S   = 86400.0_dp             ! seconds per day
  real(dp), parameter :: AU_KM   = 149597870.700_dp       ! AU in km
  real(dp), parameter :: UNIX_JD_EPOCH = 2440587.5_dp     ! JD of Unix epoch (1970-01-01)

  ! J2000 obliquity for equatorial → ecliptic rotation
  real(dp), parameter :: OBLIQUITY_DEG = 23.4392794_dp
  real(dp), parameter :: OBLIQUITY_RAD = OBLIQUITY_DEG * PI / 180.0_dp
end module constants_mod


! ─────────────────────────────────────────────────────────────────────
!  SPK reader — DAF binary format + Chebyshev evaluation
!  (from generate_observations.f90)
! ─────────────────────────────────────────────────────────────────────
module spk_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_SEGMENTS = 32

  type :: spk_segment
    real(dp) :: start_second, end_second
    integer  :: target, center, frame, data_type
    integer  :: start_i, end_i
    real(dp) :: start_jd, end_jd
    logical  :: loaded = .false.
    real(dp) :: init_epoch
    real(dp) :: intlen
    integer  :: n_intervals
    integer  :: coefficient_count
    integer  :: component_count
    real(dp), allocatable :: coeffs(:,:,:)
  end type

  type :: spk_kernel
    integer :: unit_num = -1
    integer :: n_segments = 0
    type(spk_segment) :: segments(MAX_SEGMENTS)
  end type

contains

  subroutine spk_open(filename, kernel)
    character(len=*), intent(in) :: filename
    type(spk_kernel), intent(out) :: kernel
    integer :: u, fward, nd, ni
    character(8) :: locidw

    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='old', action='read')
    kernel%unit_num = u

    read(u) locidw
    read(u) nd
    read(u) ni
    read(u, pos=77) fward

    call parse_summaries(u, fward, nd, ni, kernel)
  end subroutine

  subroutine parse_summaries(u, fward, nd, ni, kernel)
    integer, intent(in) :: u, fward, nd, ni
    type(spk_kernel), intent(inout) :: kernel
    integer :: record_number, n_summaries, i, seg_idx
    real(dp) :: next_rec, prev_rec, nsumm_d
    real(dp) :: start_sec, end_sec
    integer :: tgt, ctr, frm, dtype, si, ei
    integer :: base_pos, summary_size, step, ctrl_size

    summary_size = nd * 8 + ni * 4
    step = summary_size
    if (mod(step, 8) /= 0) step = step + (8 - mod(step, 8))
    ctrl_size = 24

    record_number = fward
    seg_idx = 0

    do while (record_number /= 0)
      base_pos = (record_number - 1) * 1024 + 1

      read(u, pos=base_pos) next_rec, prev_rec, nsumm_d
      n_summaries = int(nsumm_d)

      do i = 0, n_summaries - 1
        seg_idx = seg_idx + 1
        if (seg_idx > MAX_SEGMENTS) then
          print *, "ERROR: too many segments"
          stop 1
        end if

        read(u, pos=base_pos + ctrl_size + i * step) start_sec, end_sec, &
             tgt, ctr, frm, dtype, si, ei

        kernel%segments(seg_idx)%start_second = start_sec
        kernel%segments(seg_idx)%end_second   = end_sec
        kernel%segments(seg_idx)%target       = tgt
        kernel%segments(seg_idx)%center       = ctr
        kernel%segments(seg_idx)%frame        = frm
        kernel%segments(seg_idx)%data_type    = dtype
        kernel%segments(seg_idx)%start_i      = si
        kernel%segments(seg_idx)%end_i        = ei
        kernel%segments(seg_idx)%start_jd     = T0 + start_sec / DAY_S
        kernel%segments(seg_idx)%end_jd       = T0 + end_sec / DAY_S
      end do

      record_number = int(next_rec)
    end do

    kernel%n_segments = seg_idx
  end subroutine

  function find_segment(kernel, center, target) result(idx)
    type(spk_kernel), intent(in) :: kernel
    integer, intent(in) :: center, target
    integer :: idx, i
    idx = -1
    do i = kernel%n_segments, 1, -1
      if (kernel%segments(i)%center == center .and. &
          kernel%segments(i)%target == target) then
        idx = i
        return
      end if
    end do
  end function

  subroutine load_segment_data(kernel, idx)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: idx
    type(spk_segment) :: seg
    real(dp) :: meta(4)
    integer :: rsize_i, n_i, coeff_count, comp_count
    integer :: u, pos, total_words
    real(dp), allocatable :: raw(:)
    integer :: rec, c, k, kk

    seg = kernel%segments(idx)
    u = kernel%unit_num

    if (seg%loaded) return

    pos = (seg%end_i - 4) * 8 + 1
    read(u, pos=pos) meta(1), meta(2), meta(3), meta(4)

    kernel%segments(idx)%init_epoch = meta(1)
    kernel%segments(idx)%intlen     = meta(2)
    rsize_i = int(meta(3))
    n_i     = int(meta(4))
    kernel%segments(idx)%n_intervals = n_i

    if (seg%data_type == 2) then
      comp_count = 3
    else
      comp_count = 6
    end if
    coeff_count = (rsize_i - 2) / comp_count
    kernel%segments(idx)%coefficient_count = coeff_count
    kernel%segments(idx)%component_count   = comp_count

    total_words = rsize_i * n_i
    allocate(raw(total_words))
    pos = (seg%start_i - 1) * 8 + 1
    read(u, pos=pos) raw

    allocate(kernel%segments(idx)%coeffs(coeff_count, comp_count, n_i))

    do rec = 1, n_i
      do c = 1, comp_count
        do k = 1, coeff_count
          kk = coeff_count - k + 1
          kernel%segments(idx)%coeffs(kk, c, rec) = &
              raw((rec-1)*rsize_i + 2 + (c-1)*coeff_count + k)
        end do
      end do
    end do

    deallocate(raw)
    kernel%segments(idx)%loaded = .true.
  end subroutine

  subroutine spk_compute(kernel, center, target, tdb_whole, tdb_frac, pos)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: center, target
    real(dp), intent(in) :: tdb_whole, tdb_frac
    real(dp), intent(out) :: pos(3)
    real(dp) :: vel(3)
    call spk_compute_and_diff(kernel, center, target, tdb_whole, tdb_frac, pos, vel)
  end subroutine

  subroutine spk_compute_and_diff(kernel, center, target, tdb_whole, tdb_frac, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: center, target
    real(dp), intent(in) :: tdb_whole, tdb_frac
    real(dp), intent(out) :: pos(3), vel(3)
    integer :: idx, n_int, cc, nc
    real(dp) :: init_e, intlen
    real(dp) :: index1, offset1, index2, offset2, index3, offset_s
    integer :: interval
    real(dp) :: s, s2
    real(dp) :: w0(3), w1(3), w2(3)
    real(dp) :: dw0(3), dw1(3), dw2(3)
    real(dp) :: wlist(100, 3)
    integer :: k

    idx = find_segment(kernel, center, target)
    if (idx < 0) then
      print *, "ERROR: segment not found for center=", center, " target=", target
      stop 1
    end if

    if (.not. kernel%segments(idx)%loaded) then
      call load_segment_data(kernel, idx)
    end if

    init_e = kernel%segments(idx)%init_epoch
    intlen = kernel%segments(idx)%intlen
    n_int  = kernel%segments(idx)%n_intervals
    cc     = kernel%segments(idx)%coefficient_count
    nc     = kernel%segments(idx)%component_count

    call divmod_dp((tdb_whole - T0) * DAY_S - init_e, intlen, index1, offset1)
    call divmod_dp(tdb_frac * DAY_S, intlen, index2, offset2)
    call divmod_dp(offset1 + offset2, intlen, index3, offset_s)
    interval = int(index1 + index2 + index3)

    if (interval == n_int) then
      interval = interval - 1
      offset_s = offset_s + intlen
    end if

    interval = interval + 1

    s = 2.0_dp * offset_s / intlen - 1.0_dp
    s2 = 2.0_dp * s

    w0 = 0.0_dp
    w1 = 0.0_dp

    do k = 1, cc - 1
      w2 = w1
      w1 = w0
      w0(1) = kernel%segments(idx)%coeffs(k, 1, interval) + s2 * w1(1) - w2(1)
      w0(2) = kernel%segments(idx)%coeffs(k, 2, interval) + s2 * w1(2) - w2(2)
      w0(3) = kernel%segments(idx)%coeffs(k, 3, interval) + s2 * w1(3) - w2(3)
      wlist(k, :) = w1
    end do

    pos(1) = kernel%segments(idx)%coeffs(cc, 1, interval) + s * w0(1) - w1(1)
    pos(2) = kernel%segments(idx)%coeffs(cc, 2, interval) + s * w0(2) - w1(2)
    pos(3) = kernel%segments(idx)%coeffs(cc, 3, interval) + s * w0(3) - w1(3)

    dw0 = 0.0_dp
    dw1 = 0.0_dp

    do k = 1, cc - 1
      dw2 = dw1
      dw1 = dw0
      dw0(1) = 2.0_dp * wlist(k, 1) + dw1(1) * s2 - dw2(1)
      dw0(2) = 2.0_dp * wlist(k, 2) + dw1(2) * s2 - dw2(2)
      dw0(3) = 2.0_dp * wlist(k, 3) + dw1(3) * s2 - dw2(3)
    end do

    vel(1) = w0(1) + s * dw0(1) - dw1(1)
    vel(2) = w0(2) + s * dw0(2) - dw1(2)
    vel(3) = w0(3) + s * dw0(3) - dw1(3)

    vel = vel / intlen * 2.0_dp * DAY_S
  end subroutine

  subroutine divmod_dp(a, b, q, r)
    real(dp), intent(in)  :: a, b
    real(dp), intent(out) :: q, r
    q = floor(a / b)
    r = a - q * b
  end subroutine

  subroutine spk_close(kernel)
    type(spk_kernel), intent(inout) :: kernel
    if (kernel%unit_num > 0) close(kernel%unit_num)
    kernel%unit_num = -1
  end subroutine
end module spk_reader_mod


! ─────────────────────────────────────────────────────────────────────
!  Main program
! ─────────────────────────────────────────────────────────────────────
program artemis_trajectory
  use constants_mod
  use spk_reader_mod
  implicit none

  ! Time conversion constants
  integer, parameter  :: LEAP_SEC = 37           ! UTC leap seconds (2026)
  real(dp), parameter :: TT_OFFSET = (37.0_dp + 32.184_dp) / DAY_S  ! UTC→TT in days

  type(spk_kernel) :: kernel

  integer :: n_pts, launch_ts
  integer :: i, u_in, u_out
  integer :: ts_i
  real(dp) :: ox, oy, oz, vel_mag, rr
  real(dp) :: jd_utc, jd_whole, jd_tt_frac, jd_tdb_frac
  real(dp) :: moon_emb(3), earth_emb(3), moon_eq(3), moon_pos(3)
  real(dp) :: earth_dist, moon_dist, met_s
  real(dp) :: dx, dy, dz
  real(dp) :: cos_eps, sin_eps

  character(len=256) :: exe_dir, bsp_file
  integer :: slen

  ! Get directory of executable for data files
  call get_command_argument(0, exe_dir)
  slen = index(exe_dir, '/', back=.true.)
  if (slen > 0) then
    exe_dir = exe_dir(1:slen)
  else
    exe_dir = './'
  end if

  ! Precompute obliquity rotation (equatorial → ecliptic)
  cos_eps = cos(OBLIQUITY_RAD)
  sin_eps = sin(OBLIQUITY_RAD)

  ! Read Orion spacecraft data
  print *, "Reading Orion data..."
  open(newunit=u_in, file=trim(exe_dir) // 'artemis_orion.dat', status='old')
  read(u_in, *) n_pts, launch_ts
  print *, "  Records:", n_pts, "  Launch TS:", launch_ts

  ! Open SPK ephemeris
  print *, "Opening de440s.bsp..."
  call spk_open(trim(exe_dir) // 'de440s.bsp', kernel)

  ! Open output file
  open(newunit=u_out, file=trim(exe_dir) // 'artemis_trajectory.dat', status='replace')
  write(u_out, '(A)') '# Artemis II trajectory (Fortran-computed)'
  write(u_out, '(A)') '# Columns: timestamp met_s orion_x orion_y orion_z ' // &
                       'moon_x moon_y moon_z earth_dist_km moon_dist_km velocity_km_s range_rate_km_s'

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

contains

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
