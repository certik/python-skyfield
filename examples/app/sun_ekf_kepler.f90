module kepler_obs_kepler_mod
  use constants_mod
  implicit none
  integer, parameter :: NS = 6   ! EKF state dimension

contains

  ! ── Solve Kepler's equation M = E - e sin(E) via Newton ──
  function solve_kepler(M_in, ecc_in) result(Ek)
    real(dp), intent(in) :: M_in, ecc_in
    real(dp) :: Ek, M, dEk
    integer :: iter
    M = mod(M_in, TAU)
    if (M < 0.0_dp) M = M + TAU
    Ek = M
    do iter = 1, 30
      dEk = (Ek - ecc_in * sin(Ek) - M) / (1.0_dp - ecc_in * cos(Ek))
      Ek = Ek - dEk
      if (abs(dEk) < 1.0d-15) exit
    end do
  end function

  ! ── Cartesian (r, v, μ) → Keplerian elements ──
  subroutine cart_to_kepler(r, v, mu, a, ecc, inc, raan, argp, M0)
    real(dp), intent(in) :: r(3), v(3), mu
    real(dp), intent(out) :: a, ecc, inc, raan, argp, M0
    real(dp) :: h(3), n(3), e_vec(3), r_mag, v_mag, h_mag, n_mag
    real(dp) :: energy, E_anom, nu, cos_nu, rdotv

    r_mag = sqrt(sum(r**2))
    v_mag = sqrt(sum(v**2))

    ! Angular momentum h = r × v
    h(1) = r(2)*v(3) - r(3)*v(2)
    h(2) = r(3)*v(1) - r(1)*v(3)
    h(3) = r(1)*v(2) - r(2)*v(1)
    h_mag = sqrt(sum(h**2))

    ! Node vector n = z_hat × h
    n(1) = -h(2)
    n(2) =  h(1)
    n(3) = 0.0_dp
    n_mag = sqrt(n(1)**2 + n(2)**2)

    ! Eccentricity vector: e = (v × h)/μ - r/|r|
    e_vec(1) = (v(2)*h(3) - v(3)*h(2)) / mu - r(1) / r_mag
    e_vec(2) = (v(3)*h(1) - v(1)*h(3)) / mu - r(2) / r_mag
    e_vec(3) = (v(1)*h(2) - v(2)*h(1)) / mu - r(3) / r_mag
    ecc = sqrt(sum(e_vec**2))

    ! Semi-major axis (vis-viva)
    energy = 0.5_dp * v_mag**2 - mu / r_mag
    a = -mu / (2.0_dp * energy)

    ! Inclination
    inc = acos(max(-1.0_dp, min(1.0_dp, h(3) / h_mag)))

    ! RAAN (Ω)
    if (n_mag > 1.0d-12) then
      raan = atan2(n(2), n(1))
    else
      raan = 0.0_dp
    end if
    if (raan < 0.0_dp) raan = raan + TAU

    ! Argument of periapsis (ω)
    if (n_mag > 1.0d-12 .and. ecc > 1.0d-12) then
      cos_nu = dot_product(n, e_vec) / (n_mag * ecc)
      cos_nu = max(-1.0_dp, min(1.0_dp, cos_nu))
      argp = acos(cos_nu)
      if (e_vec(3) < 0.0_dp) argp = TAU - argp
    else
      argp = 0.0_dp
    end if

    ! True anomaly (ν)
    rdotv = dot_product(r, v)
    if (ecc > 1.0d-12) then
      cos_nu = dot_product(e_vec, r) / (ecc * r_mag)
      cos_nu = max(-1.0_dp, min(1.0_dp, cos_nu))
      nu = acos(cos_nu)
      if (rdotv < 0.0_dp) nu = TAU - nu
    else
      nu = 0.0_dp
    end if

    ! Eccentric anomaly from true anomaly
    E_anom = atan2(sqrt(1.0_dp - ecc**2) * sin(nu), ecc + cos(nu))

    ! Mean anomaly
    M0 = E_anom - ecc * sin(E_anom)
    if (M0 < 0.0_dp) M0 = M0 + TAU
  end subroutine

  ! ── Keplerian elements → Cartesian position & velocity (J2000 eq) ──
  subroutine kepler_to_cart(a, ecc, inc, raan, argp, M, mu, r_eq, v_eq)
    real(dp), intent(in) :: a, ecc, inc, raan, argp, M, mu
    real(dp), intent(out) :: r_eq(3), v_eq(3)
    real(dp) :: E_anom, nu, r_mag, h_mag
    real(dp) :: r_pf(2), v_pf(2)
    real(dp) :: cO, sO, ci, si, cw, sw
    real(dp) :: R11, R12, R21, R22, R31, R32

    E_anom = solve_kepler(M, ecc)
    nu = atan2(sqrt(1.0_dp - ecc**2) * sin(E_anom), cos(E_anom) - ecc)
    r_mag = a * (1.0_dp - ecc * cos(E_anom))
    h_mag = sqrt(mu * a * (1.0_dp - ecc**2))

    ! Perifocal frame
    r_pf(1) = r_mag * cos(nu)
    r_pf(2) = r_mag * sin(nu)
    v_pf(1) = -(mu / h_mag) * sin(nu)
    v_pf(2) =  (mu / h_mag) * (ecc + cos(nu))

    ! Rotation matrix: perifocal → equatorial
    cO = cos(raan); sO = sin(raan)
    ci = cos(inc);   si = sin(inc)
    cw = cos(argp); sw = sin(argp)

    R11 = cO*cw - sO*sw*ci
    R12 = -cO*sw - sO*cw*ci
    R21 = sO*cw + cO*sw*ci
    R22 = -sO*sw + cO*cw*ci
    R31 = sw*si
    R32 = cw*si

    r_eq(1) = R11*r_pf(1) + R12*r_pf(2)
    r_eq(2) = R21*r_pf(1) + R22*r_pf(2)
    r_eq(3) = R31*r_pf(1) + R32*r_pf(2)

    v_eq(1) = R11*v_pf(1) + R12*v_pf(2)
    v_eq(2) = R21*v_pf(1) + R22*v_pf(2)
    v_eq(3) = R31*v_pf(1) + R32*v_pf(2)
  end subroutine

  ! ──────────── Coordinate transformation routines ────────────

  function earth_rotation_angle(jd_whole, jd_frac) result(theta)
    real(dp), intent(in) :: jd_whole, jd_frac
    real(dp) :: theta, th
    th = 0.7790572732640_dp + 0.00273781191135448_dp * &
         (jd_whole - T0 + jd_frac)
    theta = mod(mod(th, 1.0_dp) + mod(jd_whole, 1.0_dp) + jd_frac, 1.0_dp)
  end function

  function gmst_hours(jd_whole, jd_frac, jd_tdb) result(gmst)
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

  function icrs_to_j2000_bias() result(B)
    real(dp) :: B(3,3)
    real(dp) :: xi0, eta0, da0, yx, zx, xy, zy, xz, yz
    xi0  = -0.0166170_dp * ASEC2RAD
    eta0 = -0.0068192_dp * ASEC2RAD
    da0  = -0.01460_dp * ASEC2RAD
    yx = -da0;   zx = xi0
    xy =  da0;   zy = eta0
    xz = -xi0;   yz = -eta0
    B(1,1) = 1.0_dp - 0.5_dp*(yx*yx + zx*zx)
    B(1,2) = xy;  B(1,3) = xz
    B(2,1) = yx
    B(2,2) = 1.0_dp - 0.5_dp*(yx*yx + zy*zy)
    B(2,3) = yz
    B(3,1) = zx;  B(3,2) = zy
    B(3,3) = 1.0_dp - 0.5_dp*(zy*zy + zx*zx)
  end function

  function compute_precession(jd_tdb) result(P)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: P(3,3), t
    real(dp) :: eps0_as, psia, omegaa, chia
    real(dp) :: eps0, sa, ca, sb, cb, sc, cc_v, sd, cd
    t = (jd_tdb - T0) / 36525.0_dp
    eps0_as = 84381.406_dp
    psia = ((((-0.0000000951_dp*t + 0.000132851_dp)*t &
              - 0.00114045_dp)*t - 1.0790069_dp)*t + 5038.481507_dp)*t
    omegaa = ((((+0.0000003337_dp*t - 0.000000467_dp)*t &
               - 0.00772503_dp)*t + 0.0512623_dp)*t - 0.025754_dp)*t + eps0_as
    chia = ((((-0.0000000560_dp*t + 0.000170663_dp)*t &
              - 0.00121197_dp)*t - 2.3814292_dp)*t + 10.556403_dp)*t
    eps0   = eps0_as * ASEC2RAD
    psia   = psia * ASEC2RAD
    omegaa = omegaa * ASEC2RAD
    chia   = chia * ASEC2RAD
    sa = sin(eps0);     ca = cos(eps0)
    sb = sin(-psia);    cb = cos(-psia)
    sc = sin(-omegaa);  cc_v = cos(-omegaa)
    sd = sin(chia);     cd = cos(chia)
    P(1,1) = cd*cb - sb*sd*cc_v
    P(1,2) = cd*sb*ca + sd*cc_v*cb*ca - sa*sd*sc
    P(1,3) = cd*sb*sa + sd*cc_v*cb*sa + ca*sd*sc
    P(2,1) = -sd*cb - sb*cd*cc_v
    P(2,2) = -sd*sb*ca + cd*cc_v*cb*ca - sa*cd*sc
    P(2,3) = -sd*sb*sa + cd*cc_v*cb*sa + ca*cd*sc
    P(3,1) = sb*sc
    P(3,2) = -sc*cb*ca - sa*cc_v
    P(3,3) = -sc*cb*sa + cc_v*ca
  end function

  function compute_PB(jd_tdb) result(M)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: M(3,3), P(3,3), B(3,3)
    B = icrs_to_j2000_bias()
    P = compute_precession(jd_tdb)
    M = matmul(P, B)
  end function

  subroutine simple_nutation(jd_tdb, dpsi, deps, mean_ob)
    real(dp), intent(in)  :: jd_tdb
    real(dp), intent(out) :: dpsi, deps, mean_ob
    real(dp) :: t, omega_node, f, d, l, lp

    t = (jd_tdb - T0) / 36525.0_dp
    mean_ob = (84381.406_dp + (-46.836769_dp + &
              (-0.0001831_dp + 0.00200340_dp * t) * t) * t) * ASEC2RAD

    ! Fundamental arguments (radians)
    omega_node = (450160.398036_dp + (-6962890.5431_dp + &
                 (7.4722_dp + 0.007702_dp * t) * t) * t) * ASEC2RAD
    omega_node = mod(omega_node, TAU)
    l = (485868.249036_dp + (1717915923.2178_dp + &
        (31.8792_dp + 0.051635_dp * t) * t) * t) * ASEC2RAD
    lp = (1287104.79305_dp + (129596581.0481_dp + &
         (-0.5532_dp + 0.000136_dp * t) * t) * t) * ASEC2RAD
    d = (1072260.70369_dp + (1602961601.2090_dp + &
        (-6.3706_dp + 0.006593_dp * t) * t) * t) * ASEC2RAD
    f = (335779.526232_dp + (1739527262.8478_dp + &
        (-12.7512_dp - 0.001037_dp * t) * t) * t) * ASEC2RAD

    dpsi = -17.2064161_dp * sin(omega_node) &
           -  1.3170907_dp * sin(2.0_dp*(f - d + omega_node)) &
           -  0.2276413_dp * sin(2.0_dp*(f + omega_node)) &
           +  0.2074554_dp * sin(2.0_dp*omega_node) &
           -  0.1491347_dp * sin(lp) &
           +  0.0130691_dp * sin(l) &
           -  0.0120323_dp * sin(2.0_dp*(f - d) + omega_node)
    dpsi = dpsi * ASEC2RAD

    deps =  9.2052331_dp * cos(omega_node) &
          + 0.5730336_dp * cos(2.0_dp*(f - d + omega_node)) &
          + 0.0978459_dp * cos(2.0_dp*(f + omega_node)) &
          - 0.0897492_dp * cos(2.0_dp*omega_node)
    deps = deps * ASEC2RAD
  end subroutine

  function compute_NPB(jd_tdb) result(M)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: M(3,3), PB(3,3), Nmat(3,3)
    real(dp) :: dpsi, deps, mean_ob, true_ob
    real(dp) :: cobm, sobm, cobt, sobt, cpsi, spsi

    PB = compute_PB(jd_tdb)
    call simple_nutation(jd_tdb, dpsi, deps, mean_ob)
    true_ob = mean_ob + deps
    cobm = cos(mean_ob);  sobm = sin(mean_ob)
    cobt = cos(true_ob);  sobt = sin(true_ob)
    cpsi = cos(dpsi);     spsi = sin(dpsi)
    Nmat(1,1) = cpsi
    Nmat(1,2) = -spsi * cobm
    Nmat(1,3) = -spsi * sobm
    Nmat(2,1) = spsi * cobt
    Nmat(2,2) = cpsi * cobm * cobt + sobm * sobt
    Nmat(2,3) = cpsi * sobm * cobt - cobm * sobt
    Nmat(3,1) = spsi * sobt
    Nmat(3,2) = cpsi * cobm * sobt - sobm * cobt
    Nmat(3,3) = cpsi * sobm * sobt + cobm * cobt
    M = matmul(Nmat, PB)
  end function

  function gast_from_gmst(gmst_h, dpsi, mean_ob, jd_tt) result(gast_h)
    real(dp), intent(in) :: gmst_h, dpsi, mean_ob, jd_tt
    real(dp) :: gast_h, t, ee_arcsec, ee_hours
    t = (jd_tt - T0) / 36525.0_dp
    ee_arcsec = dpsi / ASEC2RAD * cos(mean_ob) &
                + 0.00264_dp * sin(125.04_dp*DEG2RAD - 1934.136_dp*DEG2RAD*t) &
                + 0.000063_dp * sin(2.0_dp*(125.04_dp*DEG2RAD - 1934.136_dp*DEG2RAD*t))
    ee_hours = ee_arcsec / 54000.0_dp
    gast_h = mod(gmst_h + ee_hours, 24.0_dp)
  end function

  ! ── Predict Sun alt/az from EKF state ──
  subroutine predict_sun_altaz(x, a_fixed, t_epoch, jd_obs, lat_deg, lon_deg, &
                                alt_deg, az_deg)
    real(dp), intent(in)  :: x(NS), a_fixed, t_epoch, jd_obs, lat_deg, lon_deg
    real(dp), intent(out) :: alt_deg, az_deg

    real(dp) :: ecc, inc, raan_v, argp_v, M0, mu
    real(dp) :: n_mean, M_now, dt_s
    real(dp) :: r_eq(3), v_eq(3), r_sun(3)
    real(dp) :: d_gcrs(3), d_date(3), d_itrs(3), d_local(3)
    real(dp) :: obs_itrs(3), obs_date(3), obs_gcrs(3)
    real(dp) :: lat, lon, sinlat, coslat, slon, clon
    real(dp) :: jd_whole, jd_frac, jd_tdb
    real(dp) :: gmst_h_v, gast_h_v, gast_rad, sg, cg
    real(dp) :: M_mat(3,3)
    real(dp) :: dpsi_n, deps_n, mean_ob_n
    real(dp) :: d_unit(3), v_obs(3), vdot, dnorm
    real(dp) :: north, east, up, tmp

    ! Unpack state
    ecc = x(1); inc = x(2); raan_v = x(3)
    argp_v = x(4); M0 = x(5); mu = x(6)

    ! Mean motion and mean anomaly at observation time
    n_mean = sqrt(mu / a_fixed**3)
    dt_s = (jd_obs - t_epoch) * DAY_S
    M_now = M0 + n_mean * dt_s

    ! Keplerian → Cartesian (Earth relative to Sun, J2000)
    call kepler_to_cart(a_fixed, ecc, inc, raan_v, argp_v, M_now, mu, r_eq, v_eq)

    ! Sun geocentric position = -(Earth relative to Sun)
    r_sun = -r_eq

    ! --- Coordinate transformations (same as kalman_sim.f90) ---
    lat = lat_deg * DEG2RAD
    lon = lon_deg * DEG2RAD
    sinlat = sin(lat); coslat = cos(lat)
    slon = sin(lon); clon = cos(lon)

    jd_whole = floor(jd_obs + 0.5_dp)
    jd_frac  = jd_obs - jd_whole
    jd_tdb = jd_obs + 69.184_dp / DAY_S

    ! GMST, NPB matrix, GAST
    gmst_h_v = gmst_hours(jd_whole, jd_frac, jd_tdb)
    M_mat = compute_NPB(jd_tdb)
    call simple_nutation(jd_tdb, dpsi_n, deps_n, mean_ob_n)
    gast_h_v = gast_from_gmst(gmst_h_v, dpsi_n, mean_ob_n, jd_tdb)
    gast_rad = gast_h_v * TAU / 24.0_dp
    sg = sin(gast_rad); cg = cos(gast_rad)

    ! Observer ITRS position
    obs_itrs(1) = R_EARTH_KM * coslat * clon
    obs_itrs(2) = R_EARTH_KM * coslat * slon
    obs_itrs(3) = R_EARTH_KM * sinlat

    ! Observer in GCRS: M^T × Rz(+GAST) × obs_itrs
    obs_date(1) =  cg*obs_itrs(1) - sg*obs_itrs(2)
    obs_date(2) =  sg*obs_itrs(1) + cg*obs_itrs(2)
    obs_date(3) = obs_itrs(3)
    obs_gcrs = matmul(transpose(M_mat), obs_date)

    ! Direction in GCRS
    d_gcrs = r_sun - obs_gcrs

    ! Aberration (first-order): Earth velocity / c
    v_obs = v_eq / C_LIGHT_KMS
    dnorm = sqrt(sum(d_gcrs**2))
    if (dnorm > 1.0d-10) then
      d_unit = d_gcrs / dnorm
      vdot = dot_product(d_unit, v_obs)
      d_gcrs = d_gcrs + dnorm * (v_obs - d_unit * vdot)
    end if

    ! GCRS → true equator of date
    d_date = matmul(M_mat, d_gcrs)

    ! True equator of date → ITRS (rotate by -GAST)
    d_itrs(1) =  cg*d_date(1) + sg*d_date(2)
    d_itrs(2) = -sg*d_date(1) + cg*d_date(2)
    d_itrs(3) = d_date(3)

    ! Rotate to local meridian
    tmp        =  clon*d_itrs(1) + slon*d_itrs(2)
    d_local(2) = -slon*d_itrs(1) + clon*d_itrs(2)
    d_local(1) = tmp
    d_local(3) = d_itrs(3)

    ! Horizon coordinates
    up    = coslat*d_local(1) + sinlat*d_local(3)
    north = -sinlat*d_local(1) + coslat*d_local(3)
    east  = d_local(2)

    dnorm = sqrt(north**2 + east**2 + up**2)
    if (dnorm < 1.0d-10) then
      alt_deg = 0.0_dp; az_deg = 0.0_dp; return
    end if

    alt_deg = asin(up / dnorm) * RAD2DEG
    az_deg  = atan2(east, north) * RAD2DEG
    if (az_deg < 0.0_dp) az_deg = az_deg + 360.0_dp
  end subroutine

end module kepler_obs_kepler_mod

! ═══════════════════════════════════════════════════════════════════════
!  Main program
! ═══════════════════════════════════════════════════════════════════════
program sun_ekf_kepler
  use constants_mod
  use spk_reader_mod
  use obs_reader_mod
  use kepler_obs_kepler_mod
  implicit none

  ! GM values (DE440)
  real(dp), parameter :: GM_SUN   = 1.32712440041279419d11
  real(dp), parameter :: GM_EARTH = 3.98600435507d5

  ! State vectors
  real(dp) :: x(NS), x_true(NS), x_init(NS)
  real(dp) :: a_fixed, t_epoch

  ! EKF covariance
  real(dp) :: P_cov(NS,NS), P_pred(NS,NS)
  real(dp) :: Q_noise(NS,NS), R_noise(2,2)

  ! SPK kernel
  type(spk_kernel) :: kernel

  ! Observations
  type(observation) :: obs_all(MAX_OBS)
  integer :: n_obs_all, n_sun
  real(dp) :: s_jd(MAX_OBS), s_alt(MAX_OBS), s_az(MAX_OBS)

  ! EKF working arrays
  real(dp) :: H(2,NS), S_mat(2,2), S_inv(2,2), K_gain(NS,2)
  real(dp) :: z_obs(2), z_pred(2), innov(2)
  real(dp) :: alt_pred, az_pred, alt_p, az_p
  real(dp) :: x_pert(NS), det_S
  real(dp) :: PHt(NS,2), KH(NS,NS), I_KH(NS,NS), I_mat(NS,NS)

  ! DE440s positions
  real(dp) :: sun_ssb(3), sun_vel(3)
  real(dp) :: emb_ssb(3), emb_vel(3)
  real(dp) :: earth_emb(3), earth_emb_vel(3)
  real(dp) :: earth_ssb(3), earth_vel(3)
  real(dp) :: r_rel(3), v_rel(3)

  ! Keplerian element computation
  real(dp) :: a_comp, ecc_comp, inc_comp, raan_comp, argp_comp, M0_comp
  real(dp) :: mu_total

  ! Verification
  real(dp) :: r_check(3), v_check(3), r_diff

  ! Loop variables
  integer :: i, j, k
  real(dp) :: lat_obs, lon_obs, jd_epoch
  real(dp) :: sigma_obs
  real(dp) :: delta(NS)
  real(dp) :: rms_alt, rms_az, win_alt, win_az
  integer :: n_proc
  real(dp) :: err_omega_deg

  ! ═══════════════════════════════════════════════════
  lat_obs = 40.0_dp
  lon_obs = 0.0_dp

  print '(A)', '════════════════════════════════════════════════════════'
  print '(A)', '  2-Body Sun-Earth EKF with Keplerian Elements'
  print '(A,F8.4,A,F8.4)', '  Observer: lat=', lat_obs, ' lon=', lon_obs
  print '(A)', '════════════════════════════════════════════════════════'

  ! ── 1. Read all observations, extract Sun-only ──
  call read_observations('observations.dat', obs_all, n_obs_all, 1.0d10)
  n_sun = 0
  do i = 1, n_obs_all
    if (obs_all(i)%body == 1) then
      n_sun = n_sun + 1
      s_jd(n_sun)  = obs_all(i)%jd
      s_alt(n_sun) = obs_all(i)%alt_obs
      s_az(n_sun)  = obs_all(i)%az_obs
    end if
  end do
  print '(A,I6,A)', '  Loaded ', n_sun, ' Sun observations'

  ! ── 2. Read DE440s → compute true Keplerian elements ──
  jd_epoch = s_jd(1)
  t_epoch  = jd_epoch

  call spk_open('de440s.bsp', kernel)
  call spk_compute_and_diff(kernel, 0, 10, real(floor(jd_epoch+0.5_dp), dp), &
       jd_epoch-real(floor(jd_epoch+0.5_dp), dp), sun_ssb, sun_vel)
  call spk_compute_and_diff(kernel, 0, 3, real(floor(jd_epoch+0.5_dp), dp), &
       jd_epoch-real(floor(jd_epoch+0.5_dp), dp), emb_ssb, emb_vel)
  call spk_compute_and_diff(kernel, 3, 399, real(floor(jd_epoch+0.5_dp), dp), &
       jd_epoch-real(floor(jd_epoch+0.5_dp), dp), earth_emb, earth_emb_vel)
  call spk_close(kernel)

  earth_ssb = emb_ssb + earth_emb
  earth_vel = emb_vel + earth_emb_vel

  ! SPK reader returns velocity in km/day; convert to km/s
  sun_vel   = sun_vel   / DAY_S
  earth_vel = earth_vel / DAY_S

  r_rel = earth_ssb - sun_ssb
  v_rel = earth_vel - sun_vel
  mu_total = GM_SUN + GM_EARTH

  call cart_to_kepler(r_rel, v_rel, mu_total, a_comp, ecc_comp, inc_comp, &
       raan_comp, argp_comp, M0_comp)

  ! Verify round-trip
  call kepler_to_cart(a_comp, ecc_comp, inc_comp, raan_comp, argp_comp, &
       M0_comp, mu_total, r_check, v_check)
  r_diff = sqrt(sum((r_check - r_rel)**2))
  print '(A,ES10.3,A)', '  Round-trip error: ', r_diff, ' km'

  a_fixed = a_comp

  ! True state vector
  x_true(1) = ecc_comp
  x_true(2) = inc_comp
  x_true(3) = raan_comp
  x_true(4) = argp_comp
  x_true(5) = M0_comp
  x_true(6) = mu_total

  print '(/,A)',          '  True Keplerian elements (J2000 equatorial):'
  print '(A,F12.8)',      '    e     = ', x_true(1)
  print '(A,F10.5,A)',    '    i     = ', x_true(2) * RAD2DEG, ' deg (obliquity)'
  print '(A,F10.5,A)',    '    Omega = ', x_true(3) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    omega = ', x_true(4) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    M0    = ', x_true(5) * RAD2DEG, ' deg'
  print '(A,ES20.12,A)',  '    mu    = ', x_true(6), ' km^3/s^2'
  print '(A,F12.1,A)',    '    a     = ', a_fixed, ' km (fixed)'

  ! ── 3. Perturb initial state ──
  x_init(1) = x_true(1) * 1.20_dp            ! e: +20%
  x_init(2) = x_true(2) + 2.0_dp * DEG2RAD   ! i: +2 deg
  x_init(3) = x_true(3) + 5.0_dp * DEG2RAD   ! Omega: +5 deg
  x_init(4) = x_true(4) + 5.0_dp * DEG2RAD   ! omega: +5 deg
  x_init(5) = x_true(5) + 3.0_dp * DEG2RAD   ! M0: +3 deg
  x_init(6) = x_true(6) * 1.005_dp           ! mu: +0.5%

  print '(/,A)', '  Perturbed initial guess:'
  print '(A,F12.8,A,F6.1,A)',  '    e     = ', x_init(1), &
       '  (', (x_init(1)/x_true(1)-1.0_dp)*100.0_dp, '%)'
  print '(A,F10.5,A,F6.2,A)',  '    i     = ', x_init(2)*RAD2DEG, &
       ' deg  (delta ', (x_init(2)-x_true(2))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)',  '    Omega = ', x_init(3)*RAD2DEG, &
       ' deg  (delta ', (x_init(3)-x_true(3))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)',  '    omega = ', x_init(4)*RAD2DEG, &
       ' deg  (delta ', (x_init(4)-x_true(4))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)',  '    M0    = ', x_init(5)*RAD2DEG, &
       ' deg  (delta ', (x_init(5)-x_true(5))*RAD2DEG, ' deg)'
  print '(A,ES20.12,A,F6.3,A)', '    mu    = ', x_init(6), &
       '  (', (x_init(6)/x_true(6)-1.0_dp)*100.0_dp, '%)'

  ! ── 4. Initialize EKF ──
  x = x_init

  P_cov = 0.0_dp
  P_cov(1,1) = (0.005_dp)**2             ! sigma_e = 0.005
  P_cov(2,2) = (3.0_dp * DEG2RAD)**2     ! sigma_i = 3 deg
  P_cov(3,3) = (8.0_dp * DEG2RAD)**2     ! sigma_Omega = 8 deg
  P_cov(4,4) = (8.0_dp * DEG2RAD)**2     ! sigma_omega = 8 deg
  P_cov(5,5) = (5.0_dp * DEG2RAD)**2     ! sigma_M0 = 5 deg
  P_cov(6,6) = (0.01_dp * mu_total)**2   ! sigma_mu = 1%

  Q_noise = 0.0_dp   ! elements are constant in 2-body

  sigma_obs = 60.0_dp / 3600.0_dp   ! 60 arcsec in degrees
  R_noise = 0.0_dp
  R_noise(1,1) = sigma_obs**2
  R_noise(2,2) = sigma_obs**2

  ! FD perturbation sizes
  delta(1) = 1.0d-7             ! e
  delta(2) = 1.0d-7             ! i (rad)
  delta(3) = 1.0d-7             ! Omega (rad)
  delta(4) = 1.0d-7             ! omega (rad)
  delta(5) = 1.0d-7             ! M0 (rad)
  delta(6) = mu_total * 1.0d-7  ! mu

  I_mat = 0.0_dp
  do i = 1, NS
    I_mat(i,i) = 1.0_dp
  end do

  ! ── 5. EKF loop ──
  print '(/,A)', '  Running EKF...'
  print '(A)', '  ──────────────────────────────────────────────────────────────────────────'
  print '(A)', '   Obs#  cumRMSa" cumRMSz"  winRMSa"  winRMSz"  e_err%   i_err(d) mu_err%'

  rms_alt = 0.0_dp
  rms_az  = 0.0_dp
  win_alt = 0.0_dp
  win_az  = 0.0_dp
  n_proc  = 0

  do k = 1, n_sun
    ! ── Predict (trivial: state unchanged, covariance grows by Q) ──
    P_pred = P_cov + Q_noise

    ! ── Predicted observation ──
    call predict_sun_altaz(x, a_fixed, t_epoch, s_jd(k), lat_obs, lon_obs, &
                            alt_pred, az_pred)
    z_obs(1)  = s_alt(k)
    z_obs(2)  = s_az(k)
    z_pred(1) = alt_pred
    z_pred(2) = az_pred

    ! Innovation with azimuth wrapping
    innov(1) = z_obs(1) - z_pred(1)
    innov(2) = z_obs(2) - z_pred(2)
    if (innov(2) >  180.0_dp) innov(2) = innov(2) - 360.0_dp
    if (innov(2) < -180.0_dp) innov(2) = innov(2) + 360.0_dp

    ! ── Jacobian by finite differences ──
    do j = 1, NS
      x_pert = x
      x_pert(j) = x_pert(j) + delta(j)
      call predict_sun_altaz(x_pert, a_fixed, t_epoch, s_jd(k), &
                              lat_obs, lon_obs, alt_p, az_p)
      H(1,j) = (alt_p - alt_pred) / delta(j)
      H(2,j) = (az_p  - az_pred)  / delta(j)
      ! Azimuth wrapping in Jacobian
      if (H(2,j) >  180.0_dp / delta(j)) H(2,j) = H(2,j) - 360.0_dp / delta(j)
      if (H(2,j) < -180.0_dp / delta(j)) H(2,j) = H(2,j) + 360.0_dp / delta(j)
    end do

    ! ── Innovation covariance S = H P H^T + R ──
    PHt = matmul(P_pred, transpose(H))
    S_mat = matmul(H, PHt) + R_noise

    ! Invert 2x2
    det_S = S_mat(1,1)*S_mat(2,2) - S_mat(1,2)*S_mat(2,1)
    if (abs(det_S) < 1.0d-30) cycle
    S_inv(1,1) =  S_mat(2,2) / det_S
    S_inv(2,2) =  S_mat(1,1) / det_S
    S_inv(1,2) = -S_mat(1,2) / det_S
    S_inv(2,1) = -S_mat(2,1) / det_S

    ! ── Kalman gain K = PHt × S^{-1} ──
    K_gain = matmul(PHt, S_inv)

    ! ── State update ──
    do j = 1, NS
      x(j) = x(j) + K_gain(j,1)*innov(1) + K_gain(j,2)*innov(2)
    end do

    ! Keep angles in [0, 2π]
    x(3) = mod(x(3), TAU); if (x(3) < 0.0_dp) x(3) = x(3) + TAU
    x(4) = mod(x(4), TAU); if (x(4) < 0.0_dp) x(4) = x(4) + TAU
    x(5) = mod(x(5), TAU); if (x(5) < 0.0_dp) x(5) = x(5) + TAU
    if (x(1) < 1.0d-8) x(1) = 1.0d-8   ! e > 0
    if (x(6) < 0.0_dp) x(6) = mu_total * 0.9_dp

    ! ── Covariance update: P = (I - K H) P_pred ──
    KH = matmul(K_gain, H)
    I_KH = I_mat - KH
    P_cov = matmul(I_KH, P_pred)
    P_cov = 0.5_dp * (P_cov + transpose(P_cov))  ! symmetrize

    ! ── Running statistics ──
    n_proc = n_proc + 1
    rms_alt = rms_alt + innov(1)**2
    rms_az  = rms_az  + innov(2)**2
    ! Windowed RMS (recent 200 obs): accumulate, report, reset
    win_alt = win_alt + innov(1)**2
    win_az  = win_az  + innov(2)**2

    if (mod(n_proc, 200) == 0 .or. n_proc == 1 .or. n_proc == n_sun) then
      j = min(n_proc, 200)
      if (n_proc == 1) j = 1
      print '(I7,2F10.1,2F10.1,F9.4,F10.5,F8.3)', &
           n_proc, &
           sqrt(rms_alt / n_proc) * 3600.0_dp, &
           sqrt(rms_az  / n_proc) * 3600.0_dp, &
           sqrt(win_alt / j) * 3600.0_dp, &
           sqrt(win_az  / j) * 3600.0_dp, &
           (x(1)/x_true(1) - 1.0_dp) * 100.0_dp, &
           (x(2) - x_true(2)) * RAD2DEG, &
           (x(6)/x_true(6) - 1.0_dp) * 100.0_dp
      win_alt = 0.0_dp; win_az = 0.0_dp
    end if
  end do

  ! ── 6. Final results ──
  print '(/,A)', '════════════════════════════════════════════════════════'
  print '(A)',   '  FINAL ESTIMATES'
  print '(A)',   '════════════════════════════════════════════════════════'
  print '(A)',   '  Parameter   Estimated         True            Error'
  print '(A)',   '  ─────────────────────────────────────────────────────'
  ! Wrap angle differences for display
  err_omega_deg = (x(3) - x_true(3)) * RAD2DEG
  if (err_omega_deg >  180.0_dp) err_omega_deg = err_omega_deg - 360.0_dp
  if (err_omega_deg < -180.0_dp) err_omega_deg = err_omega_deg + 360.0_dp

  print '(A,F12.8,4X,F12.8,4X,F8.4,A)', &
       '  e       ', x(1), x_true(1), (x(1)/x_true(1)-1.0_dp)*100.0_dp, '%'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  i       ', x(2)*RAD2DEG, ' deg', x_true(2)*RAD2DEG, ' deg', &
       (x(2)-x_true(2))*RAD2DEG, ' deg'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  Omega   ', x(3)*RAD2DEG, ' deg', x_true(3)*RAD2DEG, ' deg', &
       err_omega_deg, ' deg'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  omega   ', x(4)*RAD2DEG, ' deg', x_true(4)*RAD2DEG, ' deg', &
       (x(4)-x_true(4))*RAD2DEG, ' deg'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  M0      ', x(5)*RAD2DEG, ' deg', x_true(5)*RAD2DEG, ' deg', &
       (x(5)-x_true(5))*RAD2DEG, ' deg'
  print '(A,ES20.12)', '  mu      ', x(6)
  print '(A,ES20.12,A,F8.4,A)', '  mu_true ', x_true(6), &
       '  err: ', (x(6)/x_true(6)-1.0_dp)*100.0_dp, '%'

  print '(/,A)', '  Derived quantities:'
  print '(A,F10.5,A,F10.5,A)', '    Obliquity  = ', x(2)*RAD2DEG, &
       ' deg  (true: ', x_true(2)*RAD2DEG, ' deg)'
  print '(A,F12.3,A)', '    Period     = ', &
       TAU / sqrt(x(6) / a_fixed**3) / DAY_S, ' days'
  print '(A,F12.3,A)', '    (true)     = ', &
       TAU / sqrt(x_true(6) / a_fixed**3) / DAY_S, ' days'

  print '(/,A)', '  1-sigma uncertainties (from covariance):'
  print '(A,ES10.3)',      '    sigma_e     = ', sqrt(max(0.0_dp, P_cov(1,1)))
  print '(A,F8.4,A)',      '    sigma_i     = ', &
       sqrt(max(0.0_dp, P_cov(2,2)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',      '    sigma_Omega = ', &
       sqrt(max(0.0_dp, P_cov(3,3)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',      '    sigma_omega = ', &
       sqrt(max(0.0_dp, P_cov(4,4)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',      '    sigma_M0    = ', &
       sqrt(max(0.0_dp, P_cov(5,5)))*RAD2DEG, ' deg'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_mu    = ', &
       sqrt(max(0.0_dp, P_cov(6,6))), &
       '  (', sqrt(max(0.0_dp, P_cov(6,6)))/x(6)*100.0_dp, '%)'

  print '(/,A)', '  Note: "True" elements are osculating at epoch. The EKF'
  print '(A)',   '  estimates best-fit MEAN elements over the full arc.'
  print '(A)',   '  The e difference (~5%) is expected: planetary perturbations'
  print '(A)',   '  make the osculating e vary over the year.'
  print '(A)', '════════════════════════════════════════════════════════'

end program sun_ekf_kepler
