! ─────────────────────────────────────────────────────────────────────
!  eop_mod — Earth Orientation Parameters
!
!  Supports two file formats:
!    • IERS finals2000A.data  — load_eop() / get_eop()
!    • JPL EOP2 (long/short)  — load_jpl_eop() / get_jpl_eop()
!
!  JPL EOP2 is the same source Horizons uses; it also provides the
!  celestial pole offsets (dX, dY) needed for sub-arcsecond accuracy.
! ─────────────────────────────────────────────────────────────────────
module eop_mod
  use constants_mod
  implicit none
  private

  integer, parameter :: MAX_EOP = 25000

  ! ── IERS finals2000A storage ──
  integer  :: n_eop = 0
  real(dp) :: eop_mjd(MAX_EOP)
  real(dp) :: eop_xp(MAX_EOP)       ! arcseconds
  real(dp) :: eop_yp(MAX_EOP)       ! arcseconds
  real(dp) :: eop_ut1_utc(MAX_EOP)  ! seconds

  ! ── JPL EOP2 storage ──
  integer  :: n_jpl = 0
  real(dp) :: jpl_mjd(MAX_EOP)
  real(dp) :: jpl_xp(MAX_EOP)        ! arcseconds
  real(dp) :: jpl_yp(MAX_EOP)        ! arcseconds
  real(dp) :: jpl_tai_ut1(MAX_EOP)   ! seconds
  real(dp) :: jpl_dX(MAX_EOP)        ! milliarcseconds
  real(dp) :: jpl_dY(MAX_EOP)        ! milliarcseconds

  public :: load_eop, get_eop
  public :: load_jpl_eop, get_jpl_eop

contains

  ! ════════════════════════════════════════════════════════════════════
  !  IERS finals2000A.data reader
  ! ════════════════════════════════════════════════════════════════════

  subroutine load_eop(filename)
    character(*), intent(in) :: filename
    character(len=200) :: line
    integer  :: iu, ios
    real(dp) :: mjd, xp, yp, dut

    open(newunit=iu, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(A,A)', 'ERROR: cannot open EOP file: ', trim(filename)
      error stop 1
    end if

    n_eop = 0
    do
      read(iu, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) < 68) cycle

      ! Parse fixed-width columns from finals2000A.data:
      !   cols  8-15 : MJD
      !   cols 19-27 : PM-x (arcsec)
      !   cols 38-46 : PM-y (arcsec)
      !   cols 59-68 : UT1-UTC (seconds)
      read(line(8:15),  *, iostat=ios) mjd;  if (ios /= 0) cycle
      read(line(19:27), *, iostat=ios) xp;   if (ios /= 0) cycle
      read(line(38:46), *, iostat=ios) yp;   if (ios /= 0) cycle
      read(line(59:68), *, iostat=ios) dut;  if (ios /= 0) cycle

      n_eop = n_eop + 1
      if (n_eop > MAX_EOP) then
        print '(A)', 'ERROR: EOP table overflow'
        error stop 1
      end if
      eop_mjd(n_eop)     = mjd
      eop_xp(n_eop)      = xp
      eop_yp(n_eop)      = yp
      eop_ut1_utc(n_eop) = dut
    end do
    close(iu)
  end subroutine

  subroutine get_eop(mjd, xp_as, yp_as, ut1_utc_s)
    real(dp), intent(in)  :: mjd
    real(dp), intent(out) :: xp_as, yp_as, ut1_utc_s
    integer  :: i
    real(dp) :: t

    if (n_eop == 0) then
      print '(A)', 'ERROR: EOP data not loaded'
      error stop 1
    end if

    ! Clamp to table range
    if (mjd <= eop_mjd(1)) then
      xp_as = eop_xp(1); yp_as = eop_yp(1); ut1_utc_s = eop_ut1_utc(1)
      return
    end if
    if (mjd >= eop_mjd(n_eop)) then
      xp_as = eop_xp(n_eop); yp_as = eop_yp(n_eop); ut1_utc_s = eop_ut1_utc(n_eop)
      return
    end if

    ! Binary search for bracketing interval
    i = binary_search(eop_mjd, n_eop, mjd)
    t = (mjd - eop_mjd(i)) / (eop_mjd(i+1) - eop_mjd(i))
    xp_as     = eop_xp(i)      + t * (eop_xp(i+1)      - eop_xp(i))
    yp_as     = eop_yp(i)      + t * (eop_yp(i+1)      - eop_yp(i))
    ut1_utc_s = eop_ut1_utc(i) + t * (eop_ut1_utc(i+1) - eop_ut1_utc(i))
  end subroutine

  ! ════════════════════════════════════════════════════════════════════
  !  JPL EOP2 reader  (latest_eop2.long / latest_eop2.short)
  !
  !  CSV format after "EOP2=" header line:
  !    MJD(TAI), PMx(mas), PMy(mas), TAI-UT1(ms), ...sigmas...,
  !    DX(mas), DY(mas), ...
  !
  !  Stores: xp/yp in arcsec, TAI-UT1 in seconds, dX/dY in mas.
  !  Caller computes delta_T = 32.184 + TAI_UT1.
  ! ════════════════════════════════════════════════════════════════════

  subroutine load_jpl_eop(filename)
    character(*), intent(in) :: filename
    character(len=400) :: line
    integer  :: iu, ios
    real(dp) :: mjd, pmx, pmy, tai_ut1, sig1, sig2, sig3
    real(dp) :: c1, c2, c3, dx, dy, dsig1, dsig2, dcorr
    logical  :: in_data

    open(newunit=iu, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(A,A)', 'ERROR: cannot open JPL EOP2 file: ', trim(filename)
      error stop 1
    end if

    n_jpl = 0
    in_data = .false.
    do
      read(iu, '(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Data begins after the "EOP2=" line
      if (.not. in_data) then
        if (line(1:5) == 'EOP2=' .or. line(2:6) == 'EOP2=') then
          in_data = .true.
        end if
        cycle
      end if

      ! Parse CSV: 15 comma-separated fields
      read(line, *, iostat=ios) mjd, pmx, pmy, tai_ut1, &
           sig1, sig2, sig3, c1, c2, c3, dx, dy, dsig1, dsig2, dcorr
      if (ios /= 0) cycle

      n_jpl = n_jpl + 1
      if (n_jpl > MAX_EOP) then
        print '(A)', 'ERROR: JPL EOP2 table overflow'
        error stop 1
      end if
      jpl_mjd(n_jpl)     = mjd
      jpl_xp(n_jpl)      = pmx / 1000.0_dp      ! mas → arcsec
      jpl_yp(n_jpl)      = pmy / 1000.0_dp      ! mas → arcsec
      jpl_tai_ut1(n_jpl) = tai_ut1 / 1000.0_dp  ! ms  → seconds
      jpl_dX(n_jpl)      = dx                    ! milliarcseconds
      jpl_dY(n_jpl)      = dy                    ! milliarcseconds
    end do
    close(iu)
  end subroutine

  subroutine get_jpl_eop(mjd, xp_as, yp_as, tai_ut1_s, dX_mas, dY_mas)
    real(dp), intent(in)  :: mjd
    real(dp), intent(out) :: xp_as, yp_as, tai_ut1_s, dX_mas, dY_mas
    integer  :: i
    real(dp) :: t

    if (n_jpl == 0) then
      print '(A)', 'ERROR: JPL EOP2 data not loaded'
      error stop 1
    end if

    ! Clamp to table range
    if (mjd <= jpl_mjd(1)) then
      xp_as = jpl_xp(1); yp_as = jpl_yp(1)
      tai_ut1_s = jpl_tai_ut1(1)
      dX_mas = jpl_dX(1); dY_mas = jpl_dY(1)
      return
    end if
    if (mjd >= jpl_mjd(n_jpl)) then
      xp_as = jpl_xp(n_jpl); yp_as = jpl_yp(n_jpl)
      tai_ut1_s = jpl_tai_ut1(n_jpl)
      dX_mas = jpl_dX(n_jpl); dY_mas = jpl_dY(n_jpl)
      return
    end if

    i = binary_search(jpl_mjd, n_jpl, mjd)
    t = (mjd - jpl_mjd(i)) / (jpl_mjd(i+1) - jpl_mjd(i))
    xp_as     = jpl_xp(i)      + t * (jpl_xp(i+1)      - jpl_xp(i))
    yp_as     = jpl_yp(i)      + t * (jpl_yp(i+1)      - jpl_yp(i))
    tai_ut1_s = jpl_tai_ut1(i) + t * (jpl_tai_ut1(i+1) - jpl_tai_ut1(i))
    dX_mas    = jpl_dX(i)      + t * (jpl_dX(i+1)      - jpl_dX(i))
    dY_mas    = jpl_dY(i)      + t * (jpl_dY(i+1)      - jpl_dY(i))
  end subroutine

  ! ════════════════════════════════════════════════════════════════════
  !  Shared utilities
  ! ════════════════════════════════════════════════════════════════════

  function binary_search(arr, n, val) result(idx)
    real(dp), intent(in) :: arr(:), val
    integer,  intent(in) :: n
    integer :: idx, lo, hi, mid
    lo = 1; hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (arr(mid) <= val) then
        lo = mid
      else
        hi = mid
      end if
    end do
    idx = lo
  end function

end module eop_mod
