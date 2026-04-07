! ─────────────────────────────────────────────────────────────────────
!  eop_mod — IERS Earth Orientation Parameters (finals2000A.data)
!
!  Provides interpolated polar motion (xp, yp) and UT1-UTC for any
!  Modified Julian Date within the file's range.
! ─────────────────────────────────────────────────────────────────────
module eop_mod
  use constants_mod
  implicit none
  private

  integer, parameter :: MAX_EOP = 20000
  integer  :: n_eop = 0
  real(dp) :: eop_mjd(MAX_EOP)
  real(dp) :: eop_xp(MAX_EOP)       ! arcseconds
  real(dp) :: eop_yp(MAX_EOP)       ! arcseconds
  real(dp) :: eop_ut1_utc(MAX_EOP)  ! seconds

  public :: load_eop, get_eop

contains

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
    i = binary_search(mjd)
    t = (mjd - eop_mjd(i)) / (eop_mjd(i+1) - eop_mjd(i))
    xp_as     = eop_xp(i)      + t * (eop_xp(i+1)      - eop_xp(i))
    yp_as     = eop_yp(i)      + t * (eop_yp(i+1)      - eop_yp(i))
    ut1_utc_s = eop_ut1_utc(i) + t * (eop_ut1_utc(i+1) - eop_ut1_utc(i))
  end subroutine

  function binary_search(mjd) result(idx)
    real(dp), intent(in) :: mjd
    integer :: idx, lo, hi, mid
    lo = 1; hi = n_eop
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (eop_mjd(mid) <= mjd) then
        lo = mid
      else
        hi = mid
      end if
    end do
    idx = lo
  end function

end module eop_mod
