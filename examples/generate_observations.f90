! ═══════════════════════════════════════════════════════════════════════
!  generate_observations.f90
!
!  Generates synthetic backyard Sun & Moon observations (alt, az) with
!  realistic noise for use as Kalman-filter input to an N-body sim.
!
!  For each day in the time span:
!    - Finds sunrise/sunset and moonrise/moonset by bisection
!    - Observes at rise, set, and every ~2 h while above horizon
!    - Adds systematic bias + random Gaussian noise (< 1 arcmin mean)
!
!  Input:  lat/lon (hardcoded or via source edit), de440s.bsp, nutation.dat
!  Output: observations.dat
!
!  Compile:  flang generate_observations.f90 -o generate_observations
!  Run:      ./generate_observations
! ═══════════════════════════════════════════════════════════════════════

! ─────────────────────────────────────────────────────────────────────
!  Module 1: Physical and mathematical constants
! ─────────────────────────────────────────────────────────────────────
module constants_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI      = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: TAU     = 2.0_dp * PI
  real(dp), parameter :: T0      = 2451545.0_dp           ! J2000.0 epoch
  real(dp), parameter :: DAY_S   = 86400.0_dp             ! seconds per day
  real(dp), parameter :: C_AUDAY = 173.14463267424034_dp  ! speed of light AU/day
  real(dp), parameter :: AU_M    = 149597870700.0_dp      ! AU in metres
  real(dp), parameter :: AU_KM   = AU_M / 1000.0_dp
  real(dp), parameter :: ERAD    = 6378136.6_dp           ! Earth equatorial radius (m)
  real(dp), parameter :: ANGVEL  = 7.2921150d-5           ! Earth rotation rate rad/s
  real(dp), parameter :: ASEC2RAD = 4.848136811095359935899141d-6
  real(dp), parameter :: ASEC360  = 1296000.0_dp
  real(dp), parameter :: GS      = 1.32712440017987d+20   ! GM_sun m^3/s^2
  real(dp), parameter :: C_SI    = 299792458.0_dp

  ! WGS84
  real(dp), parameter :: WGS84_RADIUS = 6378137.0_dp
  real(dp), parameter :: WGS84_INVF   = 298.257223563_dp

  ! Reciprocal masses for deflection
  real(dp), parameter :: RMASS_SUN     = 1.0_dp
  real(dp), parameter :: RMASS_JUPITER = 1047.3486_dp
  real(dp), parameter :: RMASS_SATURN  = 3497.898_dp
  real(dp), parameter :: RMASS_EARTH   = 332946.050895_dp

  real(dp), parameter :: TENTH_USEC_2_RAD = ASEC2RAD / 1.0d7

  ! Body radii
  real(dp), parameter :: SOLAR_RADIUS_KM = 695700.0_dp
  real(dp), parameter :: MOON_RADIUS_KM  = 1737.4_dp
end module constants_mod


! ─────────────────────────────────────────────────────────────────────
!  Module 2: Small linear-algebra helpers (3×3 matrices, 3-vectors)
! ─────────────────────────────────────────────────────────────────────
module linalg_mod
  use constants_mod
  implicit none
contains

  ! ── 3×3 matrix multiply ──
  function mat33_mul(A, B) result(C)
    real(dp), intent(in) :: A(3,3), B(3,3)
    real(dp) :: C(3,3)
    integer :: i, j, k
    C = 0.0_dp
    do j = 1, 3
      do k = 1, 3
        do i = 1, 3
          C(i,j) = C(i,j) + A(i,k) * B(k,j)
        end do
      end do
    end do
  end function

  ! ── 3×3 × 3-vector ──
  function mat33_vec(A, x) result(y)
    real(dp), intent(in) :: A(3,3), x(3)
    real(dp) :: y(3)
    integer :: i, k
    y = 0.0_dp
    do k = 1, 3
      do i = 1, 3
        y(i) = y(i) + A(i,k) * x(k)
      end do
    end do
  end function

  ! ── transpose ──
  function mat33_T(A) result(AT)
    real(dp), intent(in) :: A(3,3)
    real(dp) :: AT(3,3)
    integer :: i, j
    do i = 1, 3
      do j = 1, 3
        AT(i,j) = A(j,i)
      end do
    end do
  end function

  ! ── 3-vector dot ──
  function dot3(a, b) result(d)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: d
    d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function

  ! ── 3-vector length ──
  function vec_len(v) result(r)
    real(dp), intent(in) :: v(3)
    real(dp) :: r
    r = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
  end function

  ! ── rotation matrices ──
  function rot_x(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [1.0_dp,  0.0_dp, 0.0_dp]
    R(2,:) = [0.0_dp,  c,      -s     ]
    R(3,:) = [0.0_dp,  s,       c     ]
  end function

  function rot_y(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [ c,      0.0_dp, s     ]
    R(2,:) = [ 0.0_dp, 1.0_dp, 0.0_dp]
    R(3,:) = [-s,      0.0_dp, c     ]
  end function

  function rot_z(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [ c, -s,      0.0_dp]
    R(2,:) = [ s,  c,      0.0_dp]
    R(3,:) = [ 0.0_dp, 0.0_dp, 1.0_dp]
  end function

  ! ── to_spherical: xyz → (r, elevation, azimuth) ──
  subroutine to_spherical(xyz, r, elev, azim)
    real(dp), intent(in)  :: xyz(3)
    real(dp), intent(out) :: r, elev, azim
    real(dp) :: eps
    eps = tiny(1.0_dp)
    r = vec_len(xyz)
    elev = asin(xyz(3) / (r + eps))
    azim = mod(atan2(xyz(2), xyz(1)), TAU)
    if (azim < 0.0_dp) azim = azim + TAU
  end subroutine
end module linalg_mod


! ─────────────────────────────────────────────────────────────────────
!  Module 3: SPK reader — DAF binary format + Chebyshev evaluation
! ─────────────────────────────────────────────────────────────────────
module spk_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_SEGMENTS = 32

  ! One SPK segment
  type :: spk_segment
    real(dp) :: start_second, end_second
    integer  :: target, center, frame, data_type
    integer  :: start_i, end_i
    real(dp) :: start_jd, end_jd
    ! Chebyshev data (loaded on first use)
    logical  :: loaded = .false.
    real(dp) :: init_epoch           ! seconds from J2000
    real(dp) :: intlen               ! interval length (seconds)
    integer  :: n_intervals
    integer  :: coefficient_count
    integer  :: component_count
    ! coefficients(coefficient_count, component_count, n_intervals)
    ! stored in REVERSED order (highest degree first) for Clenshaw
    real(dp), allocatable :: coeffs(:,:,:)
  end type

  ! The SPK kernel
  type :: spk_kernel
    integer :: unit_num = -1
    integer :: n_segments = 0
    type(spk_segment) :: segments(MAX_SEGMENTS)
  end type

contains

  ! ── Open an SPK file ──
  subroutine spk_open(filename, kernel)
    character(len=*), intent(in) :: filename
    type(spk_kernel), intent(out) :: kernel
    integer :: u, fward, nd, ni
    character(8) :: locidw
    integer :: dummy_i

    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='old', action='read')
    kernel%unit_num = u

    ! Read file record
    read(u) locidw
    read(u) nd          ! should be 2
    read(u) ni          ! should be 6

    ! Skip locifn (60 bytes)
    read(u, pos=77) fward  ! bytes 76-79 (1-indexed: pos 77)

    ! Parse summary records
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
    integer :: name_record_base

    ! Each summary: nd doubles + ni ints = 2*8 + 6*4 = 40 bytes
    summary_size = nd * 8 + ni * 4
    step = summary_size
    if (mod(step, 8) /= 0) step = step + (8 - mod(step, 8))
    ctrl_size = 24  ! 3 doubles

    record_number = fward
    seg_idx = 0

    do while (record_number /= 0)
      base_pos = (record_number - 1) * 1024 + 1

      ! Read control area (3 doubles)
      read(u, pos=base_pos) next_rec, prev_rec, nsumm_d
      n_summaries = int(nsumm_d)

      do i = 0, n_summaries - 1
        seg_idx = seg_idx + 1
        if (seg_idx > MAX_SEGMENTS) then
          print *, "ERROR: too many segments"
          stop 1
        end if

        ! Read this summary
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

  ! ── Find segment by (center, target) ──
  function find_segment(kernel, center, target) result(idx)
    type(spk_kernel), intent(in) :: kernel
    integer, intent(in) :: center, target
    integer :: idx, i
    idx = -1
    ! Search from end (last matching segment, like jplephem)
    do i = kernel%n_segments, 1, -1
      if (kernel%segments(i)%center == center .and. &
          kernel%segments(i)%target == target) then
        idx = i
        return
      end if
    end do
  end function

  ! ── Load Chebyshev coefficients for a segment ──
  subroutine load_segment_data(kernel, idx)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: idx
    type(spk_segment) :: seg
    real(dp) :: meta(4)
    integer :: rsize_i, n_i, coeff_count, comp_count
    integer :: u, pos, total_words
    real(dp), allocatable :: raw(:), temp(:,:,:)
    integer :: rec, c, k, kk

    seg = kernel%segments(idx)
    u = kernel%unit_num

    if (seg%loaded) return

    ! Read metadata: last 4 words of segment
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

    ! Read all data words
    total_words = rsize_i * n_i
    allocate(raw(total_words))
    pos = (seg%start_i - 1) * 8 + 1
    read(u, pos=pos) raw

    ! Reshape: (n_i, rsize_i) → skip MID and RADIUS → (n_i, comp, coeff)
    ! Then reorder to (coeff, comp, n_i) reversed for Clenshaw
    allocate(kernel%segments(idx)%coeffs(coeff_count, comp_count, n_i))

    do rec = 1, n_i
      do c = 1, comp_count
        do k = 1, coeff_count
          ! In raw: row rec has rsize_i values
          ! After skipping 2 (MID, RADIUS): offset = (rec-1)*rsize_i + 2 + (c-1)*coeff_count + (k-1)
          ! We store reversed: kk = coeff_count - k + 1
          kk = coeff_count - k + 1
          kernel%segments(idx)%coeffs(kk, c, rec) = &
              raw((rec-1)*rsize_i + 2 + (c-1)*coeff_count + k)
        end do
      end do
    end do

    deallocate(raw)
    kernel%segments(idx)%loaded = .true.
  end subroutine

  ! ── Chebyshev evaluation (Clenshaw) — compute position only ──
  subroutine spk_compute(kernel, center, target, tdb_whole, tdb_frac, pos)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: center, target
    real(dp), intent(in) :: tdb_whole, tdb_frac
    real(dp), intent(out) :: pos(3)
    real(dp) :: vel(3)
    call spk_compute_and_diff(kernel, center, target, tdb_whole, tdb_frac, pos, vel)
  end subroutine

  ! ── Chebyshev evaluation — position and velocity ──
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
    real(dp) :: wlist(100, 3)   ! max 100 Chebyshev coefficients
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

    ! Two-part JD → interval index + offset (matching jplephem exactly)
    call divmod_dp((tdb_whole - T0) * DAY_S - init_e, intlen, index1, offset1)
    call divmod_dp(tdb_frac * DAY_S, intlen, index2, offset2)
    call divmod_dp(offset1 + offset2, intlen, index3, offset_s)
    interval = int(index1 + index2 + index3)

    ! Endpoint wrap
    if (interval == n_int) then
      interval = interval - 1
      offset_s = offset_s + intlen
    end if

    ! Convert to 1-indexed
    interval = interval + 1

    ! Normalized time
    s = 2.0_dp * offset_s / intlen - 1.0_dp
    s2 = 2.0_dp * s

    ! Clenshaw recurrence (coefficients are stored reversed: highest degree first)
    ! coeffs(1, :, interval) = highest degree
    ! coeffs(cc, :, interval) = constant term (degree 0)
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

    ! Position
    pos(1) = kernel%segments(idx)%coeffs(cc, 1, interval) + s * w0(1) - w1(1)
    pos(2) = kernel%segments(idx)%coeffs(cc, 2, interval) + s * w0(2) - w1(2)
    pos(3) = kernel%segments(idx)%coeffs(cc, 3, interval) + s * w0(3) - w1(3)

    ! Differentiation
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

    ! Convert rates to km/day
    vel = vel / intlen * 2.0_dp * DAY_S
  end subroutine

  ! ── divmod for real(dp) (matching Python divmod exactly) ──
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
!  Module 4: Nutation data and IAU 2000A computation
! ─────────────────────────────────────────────────────────────────────
module nutation_mod
  use constants_mod
  implicit none

  ! Nutation coefficient tables
  integer, allocatable  :: nals_t(:,:)      ! (678, 5)
  real(dp), allocatable :: lunsol_lon(:,:)  ! (678, 3)
  real(dp), allocatable :: lunsol_obl(:,:)  ! (678, 3)
  integer, allocatable  :: napl_t(:,:)      ! (687, 14)
  real(dp), allocatable :: nut_lon(:,:)     ! (687, 2)
  real(dp), allocatable :: nut_obl(:,:)     ! (687, 2)
  integer, allocatable  :: ke0_t(:,:)       ! (33, 14)
  integer, allocatable  :: ke1(:)           ! (14)
  real(dp), allocatable :: se0_t_0(:)       ! (33)
  real(dp), allocatable :: se0_t_1(:)       ! (33)

  ! Hardcoded from Skyfield
  real(dp), parameter :: se1_0 = -0.87d-6
  real(dp), parameter :: se1_1 = +0.00d-6

  ! Fundamental argument polynomial coefficients (5 arguments, 5 polynomial terms each)
  ! fa(arg, term): fa0..fa4 for each of the 5 fundamental arguments
  real(dp), parameter :: fa_coeff(5,5) = reshape([ &
    485868.249036_dp, 1287104.79305_dp, 335779.526232_dp, 1072260.70369_dp, 450160.398036_dp, &
    1717915923.2178_dp, 129596581.0481_dp, 1739527262.8478_dp, 1602961601.2090_dp, -6962890.5431_dp, &
    31.8792_dp, -0.5532_dp, -12.7512_dp, -6.3706_dp, 7.4722_dp, &
    0.051635_dp, 0.000136_dp, -0.001037_dp, 0.006593_dp, 0.007702_dp, &
    -0.00024470_dp, -0.00001149_dp, 0.00000417_dp, -0.00003169_dp, -0.00005939_dp &
  ], [5, 5])

  ! Planetary anomaly constants and coefficients (14 pairs)
  real(dp), parameter :: anom_const(14) = [ &
    2.35555598_dp, 6.24006013_dp, 1.627905234_dp, 5.198466741_dp, &
    2.18243920_dp, 4.402608842_dp, 3.176146697_dp, 1.753470314_dp, &
    6.203480913_dp, 0.599546497_dp, 0.874016757_dp, 5.481293871_dp, &
    5.321159000_dp, 0.02438175_dp ]

  real(dp), parameter :: anom_coeff(14) = [ &
    8328.6914269554_dp, 628.301955_dp, 8433.466158131_dp, 7771.3771468121_dp, &
    -33.757045_dp, 2608.7903141574_dp, 1021.3285546211_dp, 628.3075849991_dp, &
    334.0612426700_dp, 52.9690962641_dp, 21.3299104960_dp, 7.4781598567_dp, &
    3.8127774000_dp, 0.00000538691_dp ]

  logical :: nutation_loaded = .false.

contains

  ! ── Load nutation data from binary file ──
  subroutine load_nutation(filename)
    character(len=*), intent(in) :: filename
    integer :: u, r, c

    if (nutation_loaded) return

    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='old', action='read')

    ! 1. nals_t (int32)
    read(u) r, c
    allocate(nals_t(r, c))
    read(u) nals_t

    ! 2. lunisolar_longitude_coefficients (float64)
    read(u) r, c
    allocate(lunsol_lon(r, c))
    read(u) lunsol_lon

    ! 3. lunisolar_obliquity_coefficients (float64)
    read(u) r, c
    allocate(lunsol_obl(r, c))
    read(u) lunsol_obl

    ! 4. napl_t (int32)
    read(u) r, c
    allocate(napl_t(r, c))
    read(u) napl_t

    ! 5. nutation_coefficients_longitude (float64)
    read(u) r, c
    allocate(nut_lon(r, c))
    read(u) nut_lon

    ! 6. nutation_coefficients_obliquity (float64)
    read(u) r, c
    allocate(nut_obl(r, c))
    read(u) nut_obl

    ! 7. ke0_t (int32)
    read(u) r, c
    allocate(ke0_t(r, c))
    read(u) ke0_t

    ! 8. ke1 (int32, 1D: second dim is 0)
    read(u) r, c
    allocate(ke1(r))
    read(u) ke1

    ! 9. se0_t_0 (float64, 1D)
    read(u) r, c
    allocate(se0_t_0(r))
    read(u) se0_t_0

    ! 10. se0_t_1 (float64, 1D)
    read(u) r, c
    allocate(se0_t_1(r))
    read(u) se0_t_1

    close(u)
    nutation_loaded = .true.
  end subroutine

  ! ── 5 fundamental arguments (radians) ──
  subroutine fundamental_arguments(t, fa_out)
    real(dp), intent(in)  :: t
    real(dp), intent(out) :: fa_out(5)
    integer :: i
    real(dp) :: a

    do i = 1, 5
      a = fa_coeff(i, 5) * t
      a = (a + fa_coeff(i, 4)) * t
      a = (a + fa_coeff(i, 3)) * t
      a = (a + fa_coeff(i, 2)) * t
      a = a + fa_coeff(i, 1)
      a = mod(a, ASEC360)
      fa_out(i) = a * ASEC2RAD
    end do
  end subroutine

  ! ── IAU 2000A nutation ──
  subroutine iau2000a(jd_tt, dpsi, deps)
    real(dp), intent(in)  :: jd_tt
    real(dp), intent(out) :: dpsi, deps
    real(dp) :: t, fa(5), arg, sarg, carg
    real(dp) :: a_plan(14), arg_plan, sarg_p, carg_p
    integer :: i, j
    integer :: n_ls, n_pl

    t = (jd_tt - T0) / 36525.0_dp
    call fundamental_arguments(t, fa)

    n_ls = size(nals_t, 1)  ! 678
    n_pl = size(napl_t, 1)  ! 687

    dpsi = 0.0_dp
    deps = 0.0_dp

    ! Luni-solar nutation
    do i = 1, n_ls
      arg = 0.0_dp
      do j = 1, 5
        arg = arg + nals_t(i, j) * fa(j)
      end do
      sarg = sin(arg)
      carg = cos(arg)

      dpsi = dpsi + sarg * lunsol_lon(i, 1)
      dpsi = dpsi + sarg * lunsol_lon(i, 2) * t
      dpsi = dpsi + carg * lunsol_lon(i, 3)

      deps = deps + carg * lunsol_obl(i, 1)
      deps = deps + carg * lunsol_obl(i, 2) * t
      deps = deps + sarg * lunsol_obl(i, 3)
    end do

    ! Planetary nutation
    do i = 1, 14
      a_plan(i) = t * anom_coeff(i) + anom_const(i)
    end do
    a_plan(14) = a_plan(14) * t  ! last term is quadratic

    do i = 1, n_pl
      arg_plan = 0.0_dp
      do j = 1, 14
        arg_plan = arg_plan + napl_t(i, j) * a_plan(j)
      end do
      sarg_p = sin(arg_plan)
      carg_p = cos(arg_plan)

      dpsi = dpsi + sarg_p * nut_lon(i, 1)
      dpsi = dpsi + carg_p * nut_lon(i, 2)

      deps = deps + sarg_p * nut_obl(i, 1)
      deps = deps + carg_p * nut_obl(i, 2)
    end do
  end subroutine

  ! ── Equation of equinoxes complementary terms ──
  function eq_equinox_complement(jd_tt) result(c_terms)
    real(dp), intent(in) :: jd_tt
    real(dp) :: c_terms
    real(dp) :: t, fa(14), a, sa, ca
    integer :: i, j

    t = (jd_tt - T0) / 36525.0_dp

    fa(1) = ((485868.249036_dp + (715923.2178_dp + (31.8792_dp + (0.051635_dp + &
              (-0.00024470_dp) * t) * t) * t) * t) * ASEC2RAD &
              + mod(1325.0_dp * t, 1.0_dp) * TAU)
    fa(2) = ((1287104.793048_dp + (1292581.0481_dp + (-0.5532_dp + (0.000136_dp + &
              (-0.00001149_dp) * t) * t) * t) * t) * ASEC2RAD &
              + mod(99.0_dp * t, 1.0_dp) * TAU)
    fa(3) = ((335779.526232_dp + (295262.8478_dp + (-12.7512_dp + (-0.001037_dp + &
              (0.00000417_dp) * t) * t) * t) * t) * ASEC2RAD &
              + mod(1342.0_dp * t, 1.0_dp) * TAU)
    fa(4) = ((1072260.703692_dp + (1105601.2090_dp + (-6.3706_dp + (0.006593_dp + &
              (-0.00003169_dp) * t) * t) * t) * t) * ASEC2RAD &
              + mod(1236.0_dp * t, 1.0_dp) * TAU)
    fa(5) = ((450160.398036_dp + (-482890.5431_dp + (7.4722_dp + (0.007702_dp + &
              (-0.00005939_dp) * t) * t) * t) * t) * ASEC2RAD &
              + mod(-5.0_dp * t, 1.0_dp) * TAU)
    fa(6) = 4.402608842_dp + 2608.7903141574_dp * t
    fa(7) = 3.176146697_dp + 1021.3285546211_dp * t
    fa(8) = 1.753470314_dp + 628.3075849991_dp * t
    fa(9) = 6.203480913_dp + 334.0612426700_dp * t
    fa(10) = 0.599546497_dp + 52.9690962641_dp * t
    fa(11) = 0.874016757_dp + 21.3299104960_dp * t
    fa(12) = 5.481293872_dp + 7.4781598567_dp * t
    fa(13) = 5.311886287_dp + 3.8133035638_dp * t
    fa(14) = (0.024381750_dp + 0.00000538691_dp * t) * t

    ! Wrap to [0, tau)
    do i = 1, 14
      fa(i) = mod(fa(i), TAU)
    end do

    ! se1 term (ke1 dot fa)
    a = 0.0_dp
    do i = 1, 14
      a = a + ke1(i) * fa(i)
    end do
    c_terms = se1_0 * sin(a) + se1_1 * cos(a)
    c_terms = c_terms * t

    ! se0 terms (ke0_t dot fa)
    do j = 1, size(ke0_t, 1)
      a = 0.0_dp
      do i = 1, 14
        a = a + ke0_t(j, i) * fa(i)
      end do
      c_terms = c_terms + se0_t_0(j) * sin(a) + se0_t_1(j) * cos(a)
    end do

    c_terms = c_terms * ASEC2RAD
  end function

  ! ── Mean obliquity (Capitaine 2003) ──
  function mean_obliquity_rad(jd_tdb) result(eps)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: eps, t
    t = (jd_tdb - T0) / 36525.0_dp
    eps = ((((-0.0000000434_dp * t &
              - 0.000000576_dp) * t &
              + 0.00200340_dp) * t &
              - 0.0001831_dp) * t &
              - 46.836769_dp) * t + 84381.406_dp
    eps = eps * ASEC2RAD
  end function
end module nutation_mod


! ─────────────────────────────────────────────────────────────────────
!  Module 5: Astronomical computations
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


! ─────────────────────────────────────────────────────────────────────
!  Main program — generate synthetic backyard observations
! ─────────────────────────────────────────────────────────────────────
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
  real(dp), parameter :: BIAS_HALF = 25.0_dp / 3600.0_dp   ! ±25 arcsec in deg
  real(dp), parameter :: SIGMA_BASE = 20.0_dp / 3600.0_dp  ! 20 arcsec
  real(dp), parameter :: SIGMA_HOR  = 15.0_dp / 3600.0_dp  ! extra near horizon
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

  character(len=256) :: exe_dir, bsp_file
  integer :: slen

  ! ── Get data-file directory ──
  call get_command_argument(0, exe_dir)
  slen = index(exe_dir, '/', back=.true.)
  if (slen > 0) then
    exe_dir = exe_dir(1:slen)
  else
    exe_dir = './'
  end if

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
  call load_nutation(trim(exe_dir) // 'nutation.dat')
  if (index(trim(bsp_file), '/') > 0) then
    call spk_open(trim(bsp_file), kernel)
  else
    call spk_open(trim(exe_dir) // trim(bsp_file), kernel)
  end if

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
  open(newunit=u_out, file=trim(exe_dir) // 'observations.dat', &
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
