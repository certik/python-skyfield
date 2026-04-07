! ═══════════════════════════════════════════════════════════════════════
!  create_de441s.f90
!
!  Reads DE441 part-1 and part-2, extracts the Chebyshev intervals that
!  cover the same date range as DE440s (JD 2396752.5 – 2506352.5, i.e.
!  1849-01-26 to 2150-01-22), and writes de441s.bsp.
!
!  The output file is a valid DAF/SPK Type 2 kernel with 14 segments
!  (same bodies and structure as DE440s) but using DE441 coefficients.
!
!  Usage:
!    fpm run --compiler flang create_de441s
!
!  Requires: de441_part-1.bsp, de441_part-2.bsp  (from download_data.sh)
!  Produces: de441s.bsp
! ═══════════════════════════════════════════════════════════════════════
program create_de441s
  use constants_mod
  implicit none

  integer, parameter :: NSEG = 14

  ! Target time range (same as DE440s)
  real(dp), parameter :: OUT_START_JD  = 2396752.5_dp
  real(dp), parameter :: OUT_END_JD    = 2506352.5_dp
  real(dp), parameter :: OUT_START_SEC = (OUT_START_JD - T0) * DAY_S
  real(dp), parameter :: OUT_END_SEC   = (OUT_END_JD   - T0) * DAY_S

  ! Segment metadata read from a source file
  type :: seg_meta
    real(dp) :: start_sec, end_sec       ! summary time bounds
    integer  :: target, center, frame, data_type
    integer  :: start_i, end_i           ! word indices (1-based)
    real(dp) :: init_epoch               ! seconds from J2000
    real(dp) :: intlen                   ! seconds
    integer  :: rsize                    ! words per interval record
    integer  :: n_intervals
  end type

  type(seg_meta) :: p1(NSEG), p2(NSEG)
  integer :: u_p1, u_p2

  integer :: i

  call read_daf_metadata('de441_part-1.bsp', u_p1, p1)
  call read_daf_metadata('de441_part-2.bsp', u_p2, p2)

  ! Verify matching structure
  do i = 1, NSEG
    if (p1(i)%target /= p2(i)%target .or. p1(i)%center /= p2(i)%center .or. &
        p1(i)%rsize /= p2(i)%rsize) then
      print '(A,I0)', 'ERROR: segment mismatch at index ', i
      error stop 1
    end if
    ! intlen may differ for trivial 1-interval segments (Mercury/Venus bary)
    if (p1(i)%n_intervals > 1 .and. abs(p1(i)%intlen - p2(i)%intlen) > 1.0_dp) then
      print '(A,I0)', 'ERROR: intlen mismatch at index ', i
      error stop 1
    end if
  end do

  call write_de441s('de441s.bsp', u_p1, p1, u_p2, p2)

  close(u_p1)
  close(u_p2)

  print '(A)', 'Wrote de441s.bsp'

contains

  ! ── Read DAF file record + segment summaries + per-segment metadata ──
  subroutine read_daf_metadata(filename, u, segs)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: u
    type(seg_meta), intent(out) :: segs(NSEG)
    character(8) :: locidw
    integer :: nd, ni, fward
    integer :: record_num, n_summ, seg_idx, i
    integer :: base_pos, step
    real(dp) :: next_rec, prev_rec, nsumm_d
    real(dp) :: meta(4)

    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='old', action='read')

    read(u) locidw
    read(u) nd, ni
    read(u, pos=77) fward

    step = nd * 8 + ni * 4   ! 40 bytes per summary
    if (mod(step, 8) /= 0) step = step + (8 - mod(step, 8))

    record_num = fward
    seg_idx = 0

    do while (record_num /= 0)
      base_pos = (record_num - 1) * 1024 + 1
      read(u, pos=base_pos) next_rec, prev_rec, nsumm_d
      n_summ = int(nsumm_d)

      do i = 0, n_summ - 1
        seg_idx = seg_idx + 1
        read(u, pos=base_pos + 24 + i * step) &
            segs(seg_idx)%start_sec, segs(seg_idx)%end_sec, &
            segs(seg_idx)%target, segs(seg_idx)%center, &
            segs(seg_idx)%frame, segs(seg_idx)%data_type, &
            segs(seg_idx)%start_i, segs(seg_idx)%end_i

        ! Segment metadata: last 4 words of segment data
        read(u, pos=(segs(seg_idx)%end_i - 4) * 8 + 1) meta
        segs(seg_idx)%init_epoch   = meta(1)
        segs(seg_idx)%intlen       = meta(2)
        segs(seg_idx)%rsize        = int(meta(3))
        segs(seg_idx)%n_intervals  = int(meta(4))
      end do

      record_num = int(next_rec)
    end do

    if (seg_idx /= NSEG) then
      print '(A,I0,A,I0)', 'ERROR: expected ', NSEG, ' segments, got ', seg_idx
      error stop 1
    end if
  end subroutine

  ! ── Write the output DAF/SPK file ──────────────────────────────────
  subroutine write_de441s(filename, u_p1, p1, u_p2, p2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: u_p1, u_p2
    type(seg_meta), intent(in) :: p1(NSEG), p2(NSEG)

    integer :: u_out
    integer :: i, first_k, last_k, n_from_p1, n_from_p2, n_total
    integer :: out_start_i(NSEG), out_end_i(NSEG)
    integer :: word_ptr, summary_rec, name_rec, free_word
    real(dp) :: cutoff_sec, out_init_epoch, out_end_sec_seg
    real(dp), allocatable :: buf(:)
    character(8)    :: locidw, locfmt
    character(60)   :: locifn
    character(603)  :: prenul
    character(28)   :: ftpstr
    character(297)  :: pstnul
    character(40)   :: seg_name
    character(1024) :: blank_record

    ! ── Phase 1: compute output layout ────────────────────────────────
    word_ptr = 257   ! data starts after file record + comment area (2 records)

    do i = 1, NSEG
      ! Cutoff: use part-1 before part-2's init_epoch, part-2 from there
      cutoff_sec = p2(i)%init_epoch

      ! Intervals from part-1: those in [OUT_START_SEC, cutoff_sec)
      if (OUT_START_SEC < cutoff_sec .and. p1(i)%n_intervals > 1) then
        first_k = nint((OUT_START_SEC - p1(i)%init_epoch) / p1(i)%intlen)
        last_k  = nint((cutoff_sec    - p1(i)%init_epoch) / p1(i)%intlen) - 1
        first_k = max(first_k, 0)
        last_k  = min(last_k, p1(i)%n_intervals - 1)
        n_from_p1 = last_k - first_k + 1
      else
        n_from_p1 = 0
      end if

      ! Intervals from part-2: those in [cutoff_sec, OUT_END_SEC]
      if (OUT_END_SEC > cutoff_sec .and. p2(i)%n_intervals > 1) then
        first_k = nint((max(OUT_START_SEC, cutoff_sec) - p2(i)%init_epoch) / p2(i)%intlen)
        last_k  = nint((OUT_END_SEC - p2(i)%init_epoch) / p2(i)%intlen) - 1
        first_k = max(first_k, 0)
        last_k  = min(last_k, p2(i)%n_intervals - 1)
        n_from_p2 = last_k - first_k + 1
      else
        n_from_p2 = 0
      end if

      ! Trivial segments (Mercury/Venus barycenter offsets): 1 interval
      if (p1(i)%n_intervals == 1) then
        n_from_p1 = 0
        n_from_p2 = 0
        n_total = 1
      else
        n_total = n_from_p1 + n_from_p2
      end if

      out_start_i(i) = word_ptr
      ! Data: n_total * rsize words + 4 metadata words
      out_end_i(i) = word_ptr + p1(i)%rsize * n_total + 3
      word_ptr = out_end_i(i) + 1

      print '(A,I4,A,I3,A,I6,A,I6,A,I6)', &
          '  seg ', p1(i)%target, '->', p1(i)%center, &
          '  from_p1=', n_from_p1, '  from_p2=', n_from_p2, '  total=', n_total
    end do

    ! ── Compute file structure ────────────────────────────────────────
    summary_rec = (word_ptr - 2) / 128 + 2   ! first record after data
    name_rec = summary_rec + 1
    free_word = name_rec * 128 + 1

    ! ── Phase 2: write output file ────────────────────────────────────
    open(newunit=u_out, file=filename, access='stream', form='unformatted', &
         status='replace', action='write')

    ! Pre-fill key records with spaces
    blank_record = repeat(char(0), 1024)
    write(u_out, pos=1) blank_record                              ! file record
    write(u_out, pos=1025) blank_record                           ! comment area
    write(u_out, pos=(summary_rec - 1) * 1024 + 1) blank_record  ! summary record
    write(u_out, pos=(name_rec - 1) * 1024 + 1) blank_record     ! name record

    ! File record
    locidw = 'NAIF/DAF'
    locfmt = 'LTL-IEEE'
    locifn = repeat(' ', 60)
    locifn(1:9) = 'DE441S   '
    prenul = repeat(char(0), 603)
    ftpstr = repeat(char(0), 28)
    ftpstr(1:18) = 'FTPSTR:TEST:ENDFTP'
    pstnul = repeat(char(0), 297)

    write(u_out, pos=1)  locidw
    write(u_out, pos=9)  2
    write(u_out, pos=13) 6
    write(u_out, pos=17) locifn
    write(u_out, pos=77) summary_rec
    write(u_out, pos=81) summary_rec
    write(u_out, pos=85) free_word
    write(u_out, pos=89) locfmt
    write(u_out, pos=97) prenul
    write(u_out, pos=700) ftpstr
    write(u_out, pos=728) pstnul

    ! ── Write segment data ────────────────────────────────────────────
    do i = 1, NSEG
      cutoff_sec = p2(i)%init_epoch

      if (p1(i)%n_intervals == 1) then
        ! Trivial segment: write zeros with correct MID/RADIUS
        out_init_epoch = OUT_START_SEC
        out_end_sec_seg = OUT_END_SEC
        n_total = 1

        allocate(buf(p1(i)%rsize))
        buf = 0.0_dp
        buf(1) = (OUT_START_SEC + OUT_END_SEC) / 2.0_dp   ! MID
        buf(2) = (OUT_END_SEC - OUT_START_SEC) / 2.0_dp    ! RADIUS
        write(u_out, pos=(out_start_i(i) - 1) * 8 + 1) buf
        deallocate(buf)
      else
        ! Copy intervals from part-1
        out_init_epoch = OUT_START_SEC
        n_total = 0

        if (OUT_START_SEC < cutoff_sec) then
          first_k = nint((OUT_START_SEC - p1(i)%init_epoch) / p1(i)%intlen)
          last_k  = nint((cutoff_sec    - p1(i)%init_epoch) / p1(i)%intlen) - 1
          first_k = max(first_k, 0)
          last_k  = min(last_k, p1(i)%n_intervals - 1)
          n_from_p1 = last_k - first_k + 1

          if (n_from_p1 > 0) then
            allocate(buf(n_from_p1 * p1(i)%rsize))
            read(u_p1, pos=(p1(i)%start_i - 1 + first_k * p1(i)%rsize) * 8 + 1) buf
            write(u_out, pos=(out_start_i(i) - 1) * 8 + 1) buf
            n_total = n_from_p1
            deallocate(buf)
          end if
        end if

        ! Copy intervals from part-2
        if (OUT_END_SEC > cutoff_sec) then
          first_k = nint((max(OUT_START_SEC, cutoff_sec) - p2(i)%init_epoch) / p2(i)%intlen)
          last_k  = nint((OUT_END_SEC - p2(i)%init_epoch) / p2(i)%intlen) - 1
          first_k = max(first_k, 0)
          last_k  = min(last_k, p2(i)%n_intervals - 1)
          n_from_p2 = last_k - first_k + 1

          if (n_from_p2 > 0) then
            allocate(buf(n_from_p2 * p2(i)%rsize))
            read(u_p2, pos=(p2(i)%start_i - 1 + first_k * p2(i)%rsize) * 8 + 1) buf
            write(u_out, pos=(out_start_i(i) - 1 + n_total * p1(i)%rsize) * 8 + 1) buf
            n_total = n_total + n_from_p2
            deallocate(buf)
          end if
        end if

        out_end_sec_seg = out_init_epoch + n_total * p1(i)%intlen
      end if

      ! Write segment metadata (last 4 words)
      if (p1(i)%n_intervals == 1) then
        ! Trivial segment: intlen = full output span
        write(u_out, pos=(out_end_i(i) - 4) * 8 + 1) &
            out_init_epoch, OUT_END_SEC - OUT_START_SEC, &
            real(p1(i)%rsize, dp), real(n_total, dp)
      else
        write(u_out, pos=(out_end_i(i) - 4) * 8 + 1) &
            out_init_epoch, p1(i)%intlen, real(p1(i)%rsize, dp), real(n_total, dp)
      end if
    end do

    ! ── Write summary record ──────────────────────────────────────────
    write(u_out, pos=(summary_rec - 1) * 1024 + 1) 0.0_dp, 0.0_dp, real(NSEG, dp)

    do i = 1, NSEG
      write(u_out, pos=(summary_rec - 1) * 1024 + 24 + (i - 1) * 40 + 1) &
          OUT_START_SEC, OUT_END_SEC, &
          p1(i)%target, p1(i)%center, p1(i)%frame, p1(i)%data_type, &
          out_start_i(i), out_end_i(i)
    end do

    ! ── Write name record ─────────────────────────────────────────────
    do i = 1, NSEG
      seg_name = repeat(' ', 40)
      write(seg_name, '(A,I0,A,I0,A)') 'DE441S ', p1(i)%center, '->', &
          p1(i)%target, ' TYPE2'
      write(u_out, pos=(name_rec - 1) * 1024 + (i - 1) * 40 + 1) seg_name
    end do

    close(u_out)
  end subroutine

end program create_de441s
