! ─────────────────────────────────────────────────────────────────────
!  nutation_mod — IAU 2000A nutation (loads binary coefficient file)
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
