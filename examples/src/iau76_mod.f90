! ─────────────────────────────────────────────────────────────────────
!  iau76_mod — IAU 1976 precession & IAU 1980 nutation
!
!  Implements the classical precession-nutation framework used by JPL
!  Horizons (IAU76/80 + IAU 1994 equation of the equinoxes):
!
!    Precession:  Lieske (1979) — IAU 1976
!    Nutation:    Wahr (1981)   — IAU 1980  (106-term series)
!    Obliquity:   IAU 1980
!    GMST:        IAU 1982
!    Eq. of eq.:  IAU 1994
!
!  Reference: ERFA eraPmat76, eraNut80, eraObl80, eraGmst82, eraEqeq94
! ─────────────────────────────────────────────────────────────────────
module iau76_mod
  use constants_mod
  use linalg_mod
  implicit none
  private

  public :: iau80_nutation, obl80, prec76_matrix, gmst82_rad, &
            eqeq94, eqeq94_dpsi, gst94_rad, gst94_corrected_rad, compute_M80

contains

  ! ── IAU 1980 mean obliquity (ERFA eraObl80) ─────────────────────
  function obl80(jd_tt) result(eps0)
    real(dp), intent(in) :: jd_tt
    real(dp) :: eps0, t
    t = (jd_tt - T0) / 36525.0_dp
    eps0 = (84381.448_dp + (-46.8150_dp + (-0.00059_dp + &
           0.001813_dp * t) * t) * t) * ASEC2RAD
  end function

  ! ── IAU 1976 precession matrix (ERFA eraPmat76) ─────────────────
  !  From J2000.0 to date.  Uses Lieske (1979) Euler angles.
  !  ERFA convention: P = R3(-z) × R2(θ) × R3(-ζ)
  !  Our rot_z/rot_y negate the angle vs ERFA's R3/R2, so:
  !    P = rot_z(z) × rot_y(-θ) × rot_z(ζ)
  function prec76_matrix(jd_tt) result(P)
    real(dp), intent(in) :: jd_tt
    real(dp) :: P(3,3)
    real(dp) :: t, tas2r, w, zeta_a, z_a, theta_a

    t = (jd_tt - T0) / 36525.0_dp
    tas2r = t * ASEC2RAD
    w = 2306.2181_dp

    zeta_a  = (w + ((0.30188_dp) + 0.017998_dp * t) * t) * tas2r
    z_a     = (w + ((1.09468_dp) + 0.018203_dp * t) * t) * tas2r
    theta_a = (2004.3109_dp + ((-0.42665_dp) - 0.041833_dp * t) * t) * tas2r

    P = mat33_mul(rot_z(z_a), mat33_mul(rot_y(-theta_a), rot_z(zeta_a)))
  end function

  ! ── IAU 1980 nutation (ERFA eraNut80, 106 terms) ────────────────
  !  Returns dpsi, deps in radians.
  subroutine iau80_nutation(jd_tt, dpsi, deps)
    real(dp), intent(in)  :: jd_tt
    real(dp), intent(out) :: dpsi, deps

    real(dp) :: t, el, elp, f, d, om, dp_sum, de_sum, arg, s, c
    integer  :: j

    ! Conversion: 0.1 milliarcsecond → radians
    real(dp), parameter :: U2R = ASEC2RAD / 1.0e4_dp

    ! 106-term series: nl,nlp,nf,nd,nom, sp,spt, ce,cet
    ! Coefficients in units of 0.1 mas (sp,spt for longitude, ce,cet for obliquity)
    integer, parameter :: N_TERMS = 106
    integer,  parameter :: NFA(5,N_TERMS) = reshape([ &
      ! l, l', F, D, Ω  — terms 1-10
       0,  0,  0,  0,  1,   0,  0,  0,  0,  2,  -2,  0,  2,  0,  1, &
       2,  0, -2,  0,  0,  -2,  0,  2,  0,  2,   1, -1,  0, -1,  0, &
       0, -2,  2, -2,  1,   2,  0, -2,  0,  1,   0,  0,  2, -2,  2, &
       0,  1,  0,  0,  0, &
      ! 11-20
       0,  1,  2, -2,  2,   0, -1,  2, -2,  2,   0,  0,  2, -2,  1, &
       2,  0,  0, -2,  0,   0,  0,  2, -2,  0,   0,  2,  0,  0,  0, &
       0,  1,  0,  0,  1,   0,  2,  2, -2,  2,   0, -1,  0,  0,  1, &
      -2,  0,  0,  2,  1, &
      ! 21-30
       0, -1,  2, -2,  1,   2,  0,  0, -2,  1,   0,  1,  2, -2,  1, &
       1,  0,  0, -1,  0,   2,  1,  0, -2,  0,   0,  0, -2,  2,  1, &
       0,  1, -2,  2,  0,   0,  1,  0,  0,  2,  -1,  0,  0,  1,  1, &
       0,  1,  2, -2,  0, &
      ! 31-40
       0,  0,  2,  0,  2,   1,  0,  0,  0,  0,   0,  0,  2,  0,  1, &
       1,  0,  2,  0,  2,   1,  0,  0, -2,  0,  -1,  0,  2,  0,  2, &
       0,  0,  0,  2,  0,   1,  0,  0,  0,  1,  -1,  0,  0,  0,  1, &
      -1,  0,  2,  2,  2, &
      ! 41-50
       1,  0,  2,  0,  1,   0,  0,  2,  2,  2,   2,  0,  0,  0,  0, &
       1,  0,  2, -2,  2,   2,  0,  2,  0,  2,   0,  0,  2,  0,  0, &
      -1,  0,  2,  0,  1,  -1,  0,  0,  2,  1,   1,  0,  0, -2,  1, &
      -1,  0,  2,  2,  1, &
      ! 51-60
       1,  1,  0, -2,  0,   0,  1,  2,  0,  2,   0, -1,  2,  0,  2, &
       1,  0,  2,  2,  2,   1,  0,  0,  2,  0,   2,  0,  2, -2,  2, &
       0,  0,  0,  2,  1,   0,  0,  2,  2,  1,   1,  0,  2, -2,  1, &
       0,  0,  0, -2,  1, &
      ! 61-70
       1, -1,  0,  0,  0,   2,  0,  2,  0,  1,   0,  1,  0, -2,  0, &
       1,  0, -2,  0,  0,   0,  0,  0,  1,  0,   1,  1,  0,  0,  0, &
       1,  0,  2,  0,  0,   1, -1,  2,  0,  2,  -1, -1,  2,  2,  2, &
      -2,  0,  0,  0,  1, &
      ! 71-80
       3,  0,  2,  0,  2,   0, -1,  2,  2,  2,   1,  1,  2,  0,  2, &
      -1,  0,  2, -2,  1,   2,  0,  0,  0,  1,   1,  0,  0,  0,  2, &
       3,  0,  0,  0,  0,   0,  0,  2,  1,  2,  -1,  0,  0,  0,  2, &
       1,  0,  0, -4,  0, &
      ! 81-90
      -2,  0,  2,  2,  2,  -1,  0,  2,  4,  2,   2,  0,  0, -4,  0, &
       1,  1,  2, -2,  2,   1,  0,  2,  2,  1,  -2,  0,  2,  4,  2, &
      -1,  0,  4,  0,  2,   1, -1,  0, -2,  0,   2,  0,  2, -2,  1, &
       2,  0,  2,  2,  2, &
      ! 91-100
       1,  0,  0,  2,  1,   0,  0,  4, -2,  2,   3,  0,  2, -2,  2, &
       1,  0,  2, -2,  0,   0,  1,  2,  0,  1,  -1, -1,  0,  2,  1, &
       0,  0, -2,  0,  1,   0,  0,  2, -1,  2,   0,  1,  0,  2,  0, &
       1,  0, -2, -2,  0, &
      ! 101-106
       0, -1,  2,  0,  1,   1,  1,  0, -2,  1,   1,  0, -2,  2,  0, &
       2,  0,  0,  2,  0,   0,  0,  2,  4,  2,   0,  1,  0,  1,  0  &
      ], shape=[5, N_TERMS])

    ! Longitude sine coefficients: sp (constant), spt (t-dependent)
    ! Units: 0.1 milliarcsecond
    real(dp), parameter :: SP(N_TERMS) = [ &
      -171996.0_dp, 2062.0_dp, 46.0_dp, 11.0_dp, -3.0_dp, -3.0_dp, &
      -2.0_dp, 1.0_dp, -13187.0_dp, 1426.0_dp, &
      -517.0_dp, 217.0_dp, 129.0_dp, 48.0_dp, -22.0_dp, 17.0_dp, &
      -15.0_dp, -16.0_dp, -12.0_dp, -6.0_dp, &
      -5.0_dp, 4.0_dp, 4.0_dp, -4.0_dp, 1.0_dp, 1.0_dp, &
      -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
      -2274.0_dp, 712.0_dp, -386.0_dp, -301.0_dp, -158.0_dp, 123.0_dp, &
      63.0_dp, 63.0_dp, -58.0_dp, -59.0_dp, &
      -51.0_dp, -38.0_dp, 29.0_dp, 29.0_dp, -31.0_dp, 26.0_dp, &
      21.0_dp, 16.0_dp, -13.0_dp, -10.0_dp, &
      -7.0_dp, 7.0_dp, -7.0_dp, -8.0_dp, 6.0_dp, 6.0_dp, &
      -6.0_dp, -7.0_dp, 6.0_dp, -5.0_dp, &
      5.0_dp, -5.0_dp, -4.0_dp, 4.0_dp, -4.0_dp, -3.0_dp, &
      3.0_dp, -3.0_dp, -3.0_dp, -2.0_dp, &
      -3.0_dp, -3.0_dp, 2.0_dp, -2.0_dp, 2.0_dp, -2.0_dp, &
      2.0_dp, 2.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, -2.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
      -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, &
      -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
      -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp ]

    real(dp), parameter :: SPT(N_TERMS) = [ &
      -174.2_dp, 0.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, -1.6_dp, -3.4_dp, &
      1.2_dp, -0.5_dp, 0.1_dp, 0.0_dp, 0.0_dp, -0.1_dp, &
      0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      -0.2_dp, 0.1_dp, -0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.1_dp, -0.1_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ]

    ! Obliquity cosine coefficients: ce (constant), cet (t-dependent)
    real(dp), parameter :: CE(N_TERMS) = [ &
      92025.0_dp, -895.0_dp, -24.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 5736.0_dp, 54.0_dp, &
      224.0_dp, -95.0_dp, -70.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      9.0_dp, 7.0_dp, 6.0_dp, 3.0_dp, &
      3.0_dp, -2.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      977.0_dp, -7.0_dp, 200.0_dp, 129.0_dp, -1.0_dp, -53.0_dp, &
      -2.0_dp, -33.0_dp, 32.0_dp, 26.0_dp, &
      27.0_dp, 16.0_dp, -1.0_dp, -12.0_dp, 13.0_dp, -1.0_dp, &
      -10.0_dp, -8.0_dp, 7.0_dp, 5.0_dp, &
      0.0_dp, -3.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, -3.0_dp, &
      3.0_dp, 3.0_dp, -3.0_dp, 3.0_dp, &
      0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
      0.0_dp, -1.0_dp, -1.0_dp, 0.0_dp, &
      -1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ]

    real(dp), parameter :: CET(N_TERMS) = [ &
      8.9_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, -3.1_dp, -0.1_dp, &
      -0.6_dp, 0.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      -0.5_dp, 0.0_dp, 0.0_dp, -0.1_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ]

    t = (jd_tt - T0) / 36525.0_dp

    ! Fundamental arguments (Delaunay) — IAU 1980 expressions
    ! l: mean anomaly of Moon
    el = mod((485866.733_dp + (715922.633_dp + (31.310_dp + &
         0.064_dp * t) * t) * t) * ASEC2RAD + &
         mod(1325.0_dp * t, 1.0_dp) * TAU, TAU)

    ! l': mean anomaly of Sun
    elp = mod((1287099.804_dp + (1292581.224_dp + (-0.577_dp - &
          0.012_dp * t) * t) * t) * ASEC2RAD + &
          mod(99.0_dp * t, 1.0_dp) * TAU, TAU)

    ! F: Moon's mean argument of latitude
    f = mod((335778.877_dp + (295263.137_dp + (-13.257_dp + &
        0.011_dp * t) * t) * t) * ASEC2RAD + &
        mod(1342.0_dp * t, 1.0_dp) * TAU, TAU)

    ! D: mean elongation of Moon from Sun
    d = mod((1072261.307_dp + (1105601.328_dp + (-6.891_dp + &
        0.019_dp * t) * t) * t) * ASEC2RAD + &
        mod(1236.0_dp * t, 1.0_dp) * TAU, TAU)

    ! Ω: longitude of ascending node of Moon's orbit
    om = mod((450160.280_dp + (-482890.539_dp + (7.455_dp + &
         0.008_dp * t) * t) * t) * ASEC2RAD + &
         mod(-5.0_dp * t, 1.0_dp) * TAU, TAU)

    ! Sum the 106 nutation terms (largest last for precision)
    dp_sum = 0.0_dp
    de_sum = 0.0_dp
    do j = N_TERMS, 1, -1
      arg = real(NFA(1,j),dp) * el  + real(NFA(2,j),dp) * elp + &
            real(NFA(3,j),dp) * f   + real(NFA(4,j),dp) * d   + &
            real(NFA(5,j),dp) * om
      s = SP(j) + SPT(j) * t
      c = CE(j) + CET(j) * t
      if (s /= 0.0_dp) dp_sum = dp_sum + s * sin(arg)
      if (c /= 0.0_dp) de_sum = de_sum + c * cos(arg)
    end do

    ! Convert from 0.1 mas to radians
    dpsi = dp_sum * U2R
    deps = de_sum * U2R
  end subroutine

  ! ── GMST IAU 1982 (ERFA eraGmst82) ──────────────────────────────
  !  Input: UT1 as (jd_whole, jd_frac) split Julian Date.
  !  Returns GMST in radians.
  function gmst82_rad(jd_ut1_whole, jd_ut1_frac) result(gmst)
    real(dp), intent(in) :: jd_ut1_whole, jd_ut1_frac
    real(dp) :: gmst
    real(dp) :: d1, d2, t, f_sec

    ! IAU 1982 GMST-UT1 polynomial coefficients (in seconds of time)
    ! A is adjusted by -DAYSEC/2 because JD starts at noon
    real(dp), parameter :: A = 24110.54841_dp - 43200.0_dp
    real(dp), parameter :: B = 8640184.812866_dp
    real(dp), parameter :: C = 0.093104_dp
    real(dp), parameter :: DD = -6.2e-6_dp

    ! Seconds-of-time to radians: 2π / 86400
    real(dp), parameter :: S2R = TAU / DAY_S

    ! Arrange so d1 < d2 for precision
    if (jd_ut1_whole < jd_ut1_whole + jd_ut1_frac) then
      d1 = jd_ut1_whole
      d2 = jd_ut1_frac
    else
      d1 = jd_ut1_frac
      d2 = jd_ut1_whole
    end if

    ! Julian centuries from J2000.0 in UT1
    t = (d1 + (d2 - T0)) / 36525.0_dp

    ! Fractional part of JD(UT1), in seconds
    f_sec = DAY_S * (mod(d1, 1.0_dp) + mod(d2, 1.0_dp))

    ! GMST in radians, normalized to [0, 2π)
    gmst = mod(S2R * ((A + (B + (C + DD * t) * t) * t) + f_sec), TAU)
    if (gmst < 0.0_dp) gmst = gmst + TAU
  end function

  ! ── Equation of the equinoxes IAU 1994 (ERFA eraEqeq94) ─────────
  !  Returns EE in radians.  Internally calls iau80_nutation.
  function eqeq94(jd_tt) result(ee)
    real(dp), intent(in) :: jd_tt
    real(dp) :: ee
    real(dp) :: dpsi, deps

    call iau80_nutation(jd_tt, dpsi, deps)
    ee = eqeq94_dpsi(jd_tt, dpsi)
  end function

  ! ── Equation of the equinoxes given pre-computed dpsi ──────────
  !  Like eqeq94 but accepts dpsi as input (for GPS-corrected nutation).
  function eqeq94_dpsi(jd_tt, dpsi) result(ee)
    real(dp), intent(in) :: jd_tt, dpsi
    real(dp) :: ee
    real(dp) :: t, om, eps0

    t = (jd_tt - T0) / 36525.0_dp

    om = mod((450160.280_dp + (-482890.539_dp + (7.455_dp + &
         0.008_dp * t) * t) * t) * ASEC2RAD + &
         mod(-5.0_dp * t, 1.0_dp) * TAU, TAU)

    eps0 = obl80(jd_tt)

    ee = dpsi * cos(eps0) + ASEC2RAD * (0.00264_dp * sin(om) + &
         0.000063_dp * sin(om + om))
  end function

  ! ── GAST IAU 1994 = GMST82 + EQEQ94 (ERFA eraGst94) ────────────
  !  Returns GAST in radians.
  function gst94_rad(jd_ut1_whole, jd_ut1_frac, jd_tt) result(gast)
    real(dp), intent(in) :: jd_ut1_whole, jd_ut1_frac, jd_tt
    real(dp) :: gast
    gast = mod(gmst82_rad(jd_ut1_whole, jd_ut1_frac) + eqeq94(jd_tt), TAU)
    if (gast < 0.0_dp) gast = gast + TAU
  end function

  ! ── GAST with corrected dpsi ────────────────────────────────────
  !  Like gst94_rad but uses a corrected dpsi in the equation of equinoxes
  !  (for GPS/VLBI-corrected nutation).
  function gst94_corrected_rad(jd_ut1_whole, jd_ut1_frac, jd_tt, dpsi) result(gast)
    real(dp), intent(in) :: jd_ut1_whole, jd_ut1_frac, jd_tt, dpsi
    real(dp) :: gast
    gast = mod(gmst82_rad(jd_ut1_whole, jd_ut1_frac) + eqeq94_dpsi(jd_tt, dpsi), TAU)
    if (gast < 0.0_dp) gast = gast + TAU
  end function

  ! ── Combined M = N × P  (IAU 76/80, no ICRS bias) ─────────────
  !  Computes the precession-nutation matrix using IAU 1976 precession
  !  and IAU 1980 nutation.  No ICRS-to-J2000 frame bias is applied
  !  (IAU76/80 was defined in the FK5 system, not ICRS).
  !
  !  Optional dpsi_corr/deps_corr add nutation corrections (radians),
  !  e.g. derived from GPS/VLBI celestial pole offsets.
  !
  !  Returns dpsi, deps (corrected if dpsi_corr/deps_corr given),
  !  and mean obliquity, all in radians.
  subroutine compute_M80(jd_tt, M, dpsi_out, deps_out, mean_ob_out, &
                          dpsi_corr, deps_corr)
    use astro_mod, only: build_nutation_matrix
    real(dp), intent(in)  :: jd_tt
    real(dp), intent(out) :: M(3,3), dpsi_out, deps_out, mean_ob_out
    real(dp), intent(in), optional :: dpsi_corr, deps_corr
    real(dp) :: P(3,3), Nmat(3,3), true_ob

    P = prec76_matrix(jd_tt)
    call iau80_nutation(jd_tt, dpsi_out, deps_out)
    if (present(dpsi_corr)) dpsi_out = dpsi_out + dpsi_corr
    if (present(deps_corr)) deps_out = deps_out + deps_corr
    mean_ob_out = obl80(jd_tt)
    true_ob = mean_ob_out + deps_out
    Nmat = build_nutation_matrix(mean_ob_out, true_ob, dpsi_out)
    M = mat33_mul(Nmat, P)
  end subroutine

end module iau76_mod
