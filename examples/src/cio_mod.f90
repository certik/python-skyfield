! ─────────────────────────────────────────────────────────────────────
!  cio_mod — CIO-based GCRS→ITRS transformation (IERS 2010)
!
!  Implements the IAU 2006/2000A CIO-based Earth-rotation pipeline:
!    [TRS] = RPOM × R3(ERA) × Q × [CRS]
!
!  where Q is the celestial-to-intermediate matrix built from the CIP
!  coordinates (X,Y) and the CIO locator s (computed from the S06
!  series), ERA is the Earth Rotation Angle, and RPOM includes the
!  TIO locator s'.
!
!  Reference: ERFA eraC2t06a, eraS06, eraC2ixys, eraSp00, eraPom00
! ─────────────────────────────────────────────────────────────────────
module cio_mod
  use constants_mod
  use linalg_mod
  use nutation_mod, only: fundamental_arguments
  implicit none
  private

  public :: compute_cio_s, build_cio_matrix, tio_locator_sp, &
            cio_polar_motion, cio_itrs_rotation, compute_npb_fw, &
            polar_wobble_longitude, cio_itrs_rotation_no_pm

contains

  ! ── CIO locator s (ERFA eraS06) ─────────────────────────────────
  !  Computes s = CIO locator, positioning the CIO on the CIP equator.
  !  Uses the IAU 2006 series (compatible with IAU 2006/2000A PN model).
  !  Inputs: jd_tt  — TT Julian Date
  !          X, Y   — CIP coordinates in GCRS (radians)
  !  Returns: s in radians
  function compute_cio_s(jd_tt, X, Y) result(s)
    real(dp), intent(in) :: jd_tt, X, Y
    real(dp) :: s

    real(dp) :: t, fa(8), a
    real(dp) :: w0, w1, w2, w3, w4, w5
    integer :: i, j

    ! ── S06 series: s + XY/2 ──
    ! Each term has 8 integer multipliers (for fa(1:8)) and (sin,cos) coefficients.
    ! Stored as: nfa(8), s_coeff, c_coeff

    ! Polynomial part (µas)
    real(dp), parameter :: sp(6) = [ &
         94.00e-6_dp,  3808.65e-6_dp,  -122.68e-6_dp, &
      -72574.11e-6_dp,    27.98e-6_dp,    15.62e-6_dp ]

    ! ── t^0 terms (33 terms) ──
    integer, parameter :: n0 = 33
    integer, parameter :: nfa0(8,33) = reshape([ &
      0,0,0,0,1,0,0,0,  0,0,0,0,2,0,0,0,  0,0,2,-2,3,0,0,0, &
      0,0,2,-2,1,0,0,0,  0,0,2,-2,2,0,0,0,  0,0,2,0,3,0,0,0, &
      0,0,2,0,1,0,0,0,  0,0,0,0,3,0,0,0,  0,1,0,0,1,0,0,0, &
      0,1,0,0,-1,0,0,0,  1,0,0,0,-1,0,0,0,  1,0,0,0,1,0,0,0, &
      0,1,2,-2,3,0,0,0,  0,1,2,-2,1,0,0,0,  0,0,4,-4,4,0,0,0, &
      0,0,1,-1,1,-8,12,0,  0,0,2,0,0,0,0,0,  0,0,2,0,2,0,0,0, &
      1,0,2,0,3,0,0,0,  1,0,2,0,1,0,0,0,  0,0,2,-2,0,0,0,0, &
      0,1,-2,2,-3,0,0,0,  0,1,-2,2,-1,0,0,0,  0,0,0,0,0,8,-13,-1, &
      0,0,0,2,0,0,0,0,  2,0,-2,0,-1,0,0,0,  0,1,2,-2,2,0,0,0, &
      1,0,0,-2,1,0,0,0,  1,0,0,-2,-1,0,0,0,  0,0,4,-2,4,0,0,0, &
      0,0,2,-2,4,0,0,0,  1,0,-2,0,-3,0,0,0,  1,0,-2,0,-1,0,0,0 &
    ], [8, 33])
    real(dp), parameter :: sc0(2,33) = reshape([ &
      -2640.73e-6_dp, 0.39e-6_dp,  -63.53e-6_dp, 0.02e-6_dp, &
      -11.75e-6_dp, -0.01e-6_dp,  -11.21e-6_dp, -0.01e-6_dp, &
      4.57e-6_dp, 0.00e-6_dp,  -2.02e-6_dp, 0.00e-6_dp, &
      -1.98e-6_dp, 0.00e-6_dp,  1.72e-6_dp, 0.00e-6_dp, &
      1.41e-6_dp, 0.01e-6_dp,  1.26e-6_dp, 0.01e-6_dp, &
      0.63e-6_dp, 0.00e-6_dp,  0.63e-6_dp, 0.00e-6_dp, &
      -0.46e-6_dp, 0.00e-6_dp,  -0.45e-6_dp, 0.00e-6_dp, &
      -0.36e-6_dp, 0.00e-6_dp,  0.24e-6_dp, 0.12e-6_dp, &
      -0.32e-6_dp, 0.00e-6_dp,  -0.28e-6_dp, 0.00e-6_dp, &
      -0.27e-6_dp, 0.00e-6_dp,  -0.26e-6_dp, 0.00e-6_dp, &
      0.21e-6_dp, 0.00e-6_dp,  -0.19e-6_dp, 0.00e-6_dp, &
      -0.18e-6_dp, 0.00e-6_dp,  0.10e-6_dp, -0.05e-6_dp, &
      -0.15e-6_dp, 0.00e-6_dp,  0.14e-6_dp, 0.00e-6_dp, &
      0.14e-6_dp, 0.00e-6_dp,  -0.14e-6_dp, 0.00e-6_dp, &
      -0.14e-6_dp, 0.00e-6_dp,  -0.13e-6_dp, 0.00e-6_dp, &
      0.11e-6_dp, 0.00e-6_dp,  -0.11e-6_dp, 0.00e-6_dp, &
      -0.11e-6_dp, 0.00e-6_dp ], [2, 33])

    ! ── t^1 terms (3 terms) ──
    integer, parameter :: n1 = 3
    integer, parameter :: nfa1(8,3) = reshape([ &
      0,0,0,0,2,0,0,0,  0,0,0,0,1,0,0,0,  0,0,2,-2,3,0,0,0 &
    ], [8, 3])
    real(dp), parameter :: sc1(2,3) = reshape([ &
      -0.07e-6_dp, 3.57e-6_dp,  1.73e-6_dp, -0.03e-6_dp, &
      0.00e-6_dp, 0.48e-6_dp ], [2, 3])

    ! ── t^2 terms (25 terms) ──
    integer, parameter :: n2 = 25
    integer, parameter :: nfa2(8,25) = reshape([ &
      0,0,0,0,1,0,0,0,  0,0,2,-2,2,0,0,0,  0,0,2,0,2,0,0,0, &
      0,0,0,0,2,0,0,0,  0,1,0,0,0,0,0,0,  1,0,0,0,0,0,0,0, &
      0,1,2,-2,2,0,0,0,  0,0,2,0,1,0,0,0,  1,0,2,0,2,0,0,0, &
      0,1,-2,2,-2,0,0,0,  1,0,0,-2,0,0,0,0,  0,0,2,-2,1,0,0,0, &
      1,0,-2,0,-2,0,0,0,  0,0,0,2,0,0,0,0,  1,0,0,0,1,0,0,0, &
      1,0,-2,-2,-2,0,0,0,  1,0,0,0,-1,0,0,0,  1,0,2,0,1,0,0,0, &
      2,0,0,-2,0,0,0,0,  2,0,-2,0,-1,0,0,0,  0,0,2,2,2,0,0,0, &
      2,0,2,0,2,0,0,0,  2,0,0,0,0,0,0,0,  1,0,2,-2,2,0,0,0, &
      0,0,2,0,0,0,0,0 &
    ], [8, 25])
    real(dp), parameter :: sc2(2,25) = reshape([ &
      743.52e-6_dp, -0.17e-6_dp,  56.91e-6_dp, 0.06e-6_dp, &
      9.84e-6_dp, -0.01e-6_dp,  -8.85e-6_dp, 0.01e-6_dp, &
      -6.38e-6_dp, -0.05e-6_dp,  -3.07e-6_dp, 0.00e-6_dp, &
      2.23e-6_dp, 0.00e-6_dp,  1.67e-6_dp, 0.00e-6_dp, &
      1.30e-6_dp, 0.00e-6_dp,  0.93e-6_dp, 0.00e-6_dp, &
      0.68e-6_dp, 0.00e-6_dp,  -0.55e-6_dp, 0.00e-6_dp, &
      0.53e-6_dp, 0.00e-6_dp,  -0.27e-6_dp, 0.00e-6_dp, &
      -0.27e-6_dp, 0.00e-6_dp,  -0.26e-6_dp, 0.00e-6_dp, &
      -0.25e-6_dp, 0.00e-6_dp,  0.22e-6_dp, 0.00e-6_dp, &
      -0.21e-6_dp, 0.00e-6_dp,  0.20e-6_dp, 0.00e-6_dp, &
      0.17e-6_dp, 0.00e-6_dp,  0.13e-6_dp, 0.00e-6_dp, &
      -0.13e-6_dp, 0.00e-6_dp,  -0.12e-6_dp, 0.00e-6_dp, &
      -0.11e-6_dp, 0.00e-6_dp ], [2, 25])

    ! ── t^3 terms (4 terms) ──
    integer, parameter :: n3 = 4
    integer, parameter :: nfa3(8,4) = reshape([ &
      0,0,0,0,1,0,0,0,  0,0,2,-2,2,0,0,0, &
      0,0,2,0,2,0,0,0,  0,0,0,0,2,0,0,0 &
    ], [8, 4])
    real(dp), parameter :: sc3(2,4) = reshape([ &
      0.30e-6_dp, -23.42e-6_dp,  -0.03e-6_dp, -1.46e-6_dp, &
      -0.01e-6_dp, -0.25e-6_dp,  0.00e-6_dp, 0.23e-6_dp ], [2, 4])

    ! ── t^4 terms (1 term) ──
    integer, parameter :: n4 = 1
    integer, parameter :: nfa4(8,1) = reshape([ &
      0,0,0,0,1,0,0,0 &
    ], [8, 1])
    real(dp), parameter :: sc4(2,1) = reshape([ &
      -0.26e-6_dp, -0.01e-6_dp ], [2, 1])

    ! ── Evaluate ──
    t = (jd_tt - T0) / 36525.0_dp

    ! 8 fundamental arguments (IERS 2003)
    call fundamental_arguments(t, fa(1:5))
    fa(6) = mod(3.176146697_dp + 1021.3285546211_dp * t, TAU)  ! Venus
    fa(7) = mod(1.753470314_dp +  628.3075849991_dp * t, TAU)  ! Earth
    fa(8) = (0.024381750_dp + 0.00000538691_dp * t) * t        ! gen. prec.

    ! Accumulate series
    w0 = sp(1); w1 = sp(2); w2 = sp(3)
    w3 = sp(4); w4 = sp(5); w5 = sp(6)

    do i = n0, 1, -1
      a = 0.0_dp
      do j = 1, 8
        a = a + real(nfa0(j,i), dp) * fa(j)
      end do
      w0 = w0 + sc0(1,i) * sin(a) + sc0(2,i) * cos(a)
    end do

    do i = n1, 1, -1
      a = 0.0_dp
      do j = 1, 8
        a = a + real(nfa1(j,i), dp) * fa(j)
      end do
      w1 = w1 + sc1(1,i) * sin(a) + sc1(2,i) * cos(a)
    end do

    do i = n2, 1, -1
      a = 0.0_dp
      do j = 1, 8
        a = a + real(nfa2(j,i), dp) * fa(j)
      end do
      w2 = w2 + sc2(1,i) * sin(a) + sc2(2,i) * cos(a)
    end do

    do i = n3, 1, -1
      a = 0.0_dp
      do j = 1, 8
        a = a + real(nfa3(j,i), dp) * fa(j)
      end do
      w3 = w3 + sc3(1,i) * sin(a) + sc3(2,i) * cos(a)
    end do

    do i = n4, 1, -1
      a = 0.0_dp
      do j = 1, 8
        a = a + real(nfa4(j,i), dp) * fa(j)
      end do
      w4 = w4 + sc4(1,i) * sin(a) + sc4(2,i) * cos(a)
    end do

    ! s + XY/2
    s = (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) * ASEC2RAD &
        - X * Y / 2.0_dp
  end function

  ! ── Build CIO matrix Q from X, Y, s (ERFA eraC2ixys) ───────────
  !  Q = Rz(E) · Ry(d) · Rz(-(E+s))
  !  where E = atan2(Y, X), d = atan(sqrt(r²/(1-r²))), r² = X²+Y².
  !  Note: uses OUR rotation convention (not ERFA's).
  !  ERFA: Rz(e) means [cos(e) sin(e); -sin(e) cos(e); ...]
  !       = our rot_z(-e)
  !  So ERFA's sequence Rz(e) · Ry(d) · Rz(-(e+s)) translates to:
  !       rot_z(-e) · rot_y(-d) · rot_z(e+s)
  function build_cio_matrix(X, Y, s) result(Q)
    real(dp), intent(in) :: X, Y, s
    real(dp) :: Q(3,3)
    real(dp) :: r2, E, d

    r2 = X * X + Y * Y
    if (r2 > 0.0_dp) then
      E = atan2(Y, X)
    else
      E = 0.0_dp
    end if
    d = atan(sqrt(r2 / (1.0_dp - r2)))

    Q = mat33_mul(rot_z(E + s), mat33_mul(rot_y(-d), rot_z(-E)))
  end function

  ! ── TIO locator s' (ERFA eraSp00) ──────────────────────────────
  !  ~47 µas/century secular drift.
  function tio_locator_sp(jd_tt) result(sp)
    real(dp), intent(in) :: jd_tt
    real(dp) :: sp, t
    t = (jd_tt - T0) / 36525.0_dp
    sp = -47.0e-6_dp * t * ASEC2RAD
  end function

  ! ── CIO polar motion matrix (ERFA eraPom00) ────────────────────
  !  RPOM = Rx(-yp) · Ry(-xp) · Rz(sp)  (ERFA convention)
  !       = rot_x(yp) · rot_y(xp) · rot_z(-sp)  (our convention)
  function cio_polar_motion(xp_as, yp_as, sp_rad) result(RPOM)
    real(dp), intent(in) :: xp_as, yp_as, sp_rad
    real(dp) :: RPOM(3,3)
    real(dp) :: xp_rad, yp_rad
    xp_rad = xp_as * ASEC2RAD
    yp_rad = yp_as * ASEC2RAD
    RPOM = mat33_mul(rot_x(yp_rad), mat33_mul(rot_y(xp_rad), rot_z(-sp_rad)))
  end function

  ! ── Full CIO-based GCRS→ITRS rotation (ERFA eraC2tcio) ─────────
  !  R_itrs = RPOM × rot_z(-ERA) × Q
  function cio_itrs_rotation(Q, era_rad, RPOM) result(R)
    real(dp), intent(in) :: Q(3,3), era_rad, RPOM(3,3)
    real(dp) :: R(3,3)
    R = mat33_mul(RPOM, mat33_mul(rot_z(-era_rad), Q))
  end function

  ! ── GCRS→ITRS without polar motion (ERA × Q only) ──────────────
  !  R = rot_z(-ERA) × Q   (CIP-based terrestrial frame, not ITRS)
  function cio_itrs_rotation_no_pm(Q, era_rad) result(R)
    real(dp), intent(in) :: Q(3,3), era_rad
    real(dp) :: R(3,3)
    R = mat33_mul(rot_z(-era_rad), Q)
  end function

  ! ── Classical polar-wobble longitude correction ──────────────────
  !  Δλ = (xp·sin(λ) + yp·cos(λ)) · tan(φ)
  !
  !  This is the classical approach used by JPL Horizons: polar motion
  !  shifts the observer's effective longitude (affecting sidereal time
  !  / hour angle) but does NOT change the geodetic latitude used for
  !  the horizon frame.  The modern RPOM matrix approach corrects both
  !  latitude and longitude; the difference is ~0.14" in altitude.
  !
  !  Returns Δλ in radians.
  function polar_wobble_longitude(lat_rad, lon_rad, xp_as, yp_as) result(dlam)
    real(dp), intent(in) :: lat_rad, lon_rad, xp_as, yp_as
    real(dp) :: dlam
    real(dp) :: xp_rad, yp_rad
    xp_rad = xp_as * ASEC2RAD
    yp_rad = yp_as * ASEC2RAD
    dlam = (xp_rad * sin(lon_rad) + yp_rad * cos(lon_rad)) * tan(lat_rad)
  end function

  ! ── NPB matrix via Fukushima-Williams angles (ERFA eraPn06a) ────
  !
  !  NPB = R1(-(epsa+deps)) · R3(-(psib+dpsi)) · R1(phib) · R3(gamb)
  !
  !  In our rotation convention (ERFA Ri(θ) = rot_i(-θ)):
  !  NPB = rot_x(epsa+deps) · rot_z(psib+dpsi) · rot_x(-phib) · rot_z(-gamb)
  !
  !  FW angles from eraPfw06 (Hilton et al. 2006):
  !    gamb, phib, psib  — bias+precession angles (arcsec polynomials in t)
  !    epsa              — mean obliquity of date (IAU 2006)
  !
  !  The nutation dpsi, deps from IAU 2000A are added to psib and epsa.
  !  Also applies the IAU 2006 J2-rate correction (Nut06a).
  !
  !  Returns: NPB matrix, d_psi, d_eps (nutation in radians), mean_ob
  subroutine compute_npb_fw(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
    use nutation_mod, only: iau2000a, mean_obliquity_rad
    real(dp), intent(in)  :: jd_tt, jd_tdb
    real(dp), intent(out) :: M(3,3), d_psi, d_eps, mean_ob
    real(dp) :: t, gamb, phib, psib, epsa
    real(dp) :: dpsi_raw, deps_raw, fj2

    t = (jd_tdb - T0) / 36525.0_dp

    ! Fukushima-Williams bias+precession angles (eraPfw06)
    gamb = (    -0.052928_dp     + &
           (    10.556378_dp     + &
           (     0.4932044_dp    + &
           (    -0.00031238_dp   + &
           (    -0.000002788_dp  + &
           (     0.0000000260_dp ) &
           * t) * t) * t) * t) * t) * ASEC2RAD

    phib = ( 84381.412819_dp     + &
           (   -46.811016_dp     + &
           (     0.0511268_dp    + &
           (     0.00053289_dp   + &
           (    -0.000000440_dp  + &
           (    -0.0000000176_dp ) &
           * t) * t) * t) * t) * t) * ASEC2RAD

    psib = (    -0.041775_dp     + &
           (  5038.481484_dp     + &
           (     1.5584175_dp    + &
           (    -0.00018522_dp   + &
           (    -0.000026452_dp  + &
           (    -0.0000000148_dp ) &
           * t) * t) * t) * t) * t) * ASEC2RAD

    ! Mean obliquity (IAU 2006, same as eraPfw06 epsa)
    mean_ob = mean_obliquity_rad(jd_tdb)
    epsa = mean_ob

    ! IAU 2000A nutation
    call iau2000a(jd_tt, dpsi_raw, deps_raw)
    d_psi = dpsi_raw * TENTH_USEC_2_RAD
    d_eps = deps_raw * TENTH_USEC_2_RAD

    ! IAU 2006 J2-rate correction (eraNut06a)
    fj2 = -2.7774e-6_dp * t
    d_psi = d_psi + d_psi * (0.4697e-6_dp + fj2)
    d_eps = d_eps + d_eps * fj2

    ! Build NPB via FW: R1(-(epsa+deps)) · R3(-(psib+dpsi)) · R1(phib) · R3(gamb)
    ! In our convention: rot_x(epsa+deps) · rot_z(psib+dpsi) · rot_x(-phib) · rot_z(-gamb)
    M = mat33_mul(rot_x(epsa + d_eps), &
        mat33_mul(rot_z(psib + d_psi), &
        mat33_mul(rot_x(-phib), rot_z(-gamb))))
  end subroutine

end module cio_mod
