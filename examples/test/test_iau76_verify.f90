! Verification of iau76_mod against ERFA reference values
program test_iau76_verify
  use constants_mod
  use iau76_mod
  implicit none
  real(dp) :: dpsi, deps, eps0, gmst, gast, M(3,3), mean_ob
  real(dp) :: jd_tt
  integer :: n_fail

  n_fail = 0
  jd_tt = 2460677.500800_dp

  ! obl80: ERFA gives 84369.742638"
  eps0 = obl80(jd_tt)
  call chk('obl80', eps0 / ASEC2RAD, 84369.742638_dp, 0.001_dp, n_fail)

  ! nut80: ERFA gives dpsi=0.320406", deps=8.553826"
  call iau80_nutation(jd_tt, dpsi, deps)
  call chk('nut80 dpsi', dpsi / ASEC2RAD, 0.320406_dp, 0.001_dp, n_fail)
  call chk('nut80 deps', deps / ASEC2RAD, 8.553826_dp, 0.001_dp, n_fail)

  ! gmst82: ERFA gives 18.7522730571332 h
  gmst = gmst82_rad(2460677.0_dp, -0.000300_dp)
  call chk('gmst82', gmst * 12.0_dp / PI, 18.7522730571332_dp, 1.0e-10_dp, n_fail)

  ! gst94: We use TT for eqeq94, ERFA uses UT1 → ~1e-6 h difference
  gast = gst94_rad(2460677.0_dp, -0.000300_dp, jd_tt)
  call chk('gst94', gast * 12.0_dp / PI, 18.7522775324379_dp, 2.0e-6_dp, n_fail)

  print '(A)', ''
  if (n_fail == 0) then
    print '(A)', 'PASS — iau76_mod verified against ERFA.'
  else
    print '(I0,A)', n_fail, ' check(s) FAILED.'
    error stop 1
  end if

contains
  subroutine chk(label, got, ref, tol, n_fail)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: got, ref, tol
    integer, intent(inout) :: n_fail
    real(dp) :: diff
    diff = abs(got - ref)
    if (diff <= tol) then
      print '(A,A,A,ES10.2)', '  PASS  ', label, '  diff=', diff
    else
      print '(A,A,A,ES10.2,A,ES10.2)', '  FAIL  ', label, '  diff=', diff, '  tol=', tol
      n_fail = n_fail + 1
    end if
  end subroutine
end program
