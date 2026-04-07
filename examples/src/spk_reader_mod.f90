! ─────────────────────────────────────────────────────────────────────
!  spk_reader_mod — DAF/SPK binary reader + Chebyshev evaluation
! ─────────────────────────────────────────────────────────────────────
module spk_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_SEGMENTS = 64

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
    real(dp), allocatable :: raw(:)
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

    ! Clenshaw recurrence (coefficients stored reversed: highest degree first)
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
