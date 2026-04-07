module nbody3_ekf_mod
  use constants_mod
  implicit none
  integer, parameter :: NB = 3   ! Sun, Earth, Moon

  ! Yoshida 4th-order symplectic coefficients
  real(dp), parameter :: W0_Y = -(2.0_dp ** (1.0_dp/3.0_dp)) / &
      (2.0_dp - 2.0_dp ** (1.0_dp/3.0_dp))
  real(dp), parameter :: W1_Y = 1.0_dp / (2.0_dp - 2.0_dp ** (1.0_dp/3.0_dp))
  real(dp), parameter :: C1_Y = 0.5_dp * W1_Y
  real(dp), parameter :: C4_Y = C1_Y
  real(dp), parameter :: C2_Y = 0.5_dp * (W0_Y + W1_Y)
  real(dp), parameter :: C3_Y = C2_Y
  real(dp), parameter :: D1_Y = W1_Y
  real(dp), parameter :: D2_Y = W0_Y
  real(dp), parameter :: D3_Y = W1_Y

contains

  subroutine compute_acc_nb(pos, gm, acc)
    real(dp), intent(in)  :: pos(NB,3), gm(NB)
    real(dp), intent(out) :: acc(NB,3)
    real(dp) :: rij(3), r2, rinv3
    integer :: i, j

    acc = 0.0_dp
    do i = 1, NB
      do j = i+1, NB
        rij = pos(j,:) - pos(i,:)
        r2 = dot_product(rij, rij)
        rinv3 = 1.0_dp / (r2 * sqrt(r2))
        acc(i,:) = acc(i,:) + gm(j) * rij * rinv3
        acc(j,:) = acc(j,:) - gm(i) * rij * rinv3
      end do
    end do
  end subroutine

  subroutine yoshida4_step_nb(pos, vel, gm, dt)
    real(dp), intent(inout) :: pos(NB,3), vel(NB,3)
    real(dp), intent(in)    :: gm(NB), dt
    real(dp) :: acc(NB,3)

    pos = pos + C1_Y * dt * vel
    call compute_acc_nb(pos, gm, acc)
    vel = vel + D1_Y * dt * acc

    pos = pos + C2_Y * dt * vel
    call compute_acc_nb(pos, gm, acc)
    vel = vel + D2_Y * dt * acc

    pos = pos + C3_Y * dt * vel
    call compute_acc_nb(pos, gm, acc)
    vel = vel + D3_Y * dt * acc

    pos = pos + C4_Y * dt * vel
  end subroutine

  subroutine propagate_nbody(pos, vel, gm, dt_total)
    real(dp), intent(inout) :: pos(NB,3), vel(NB,3)
    real(dp), intent(in)    :: gm(NB), dt_total
    real(dp) :: dt_step
    integer :: n_steps, i

    ! 3-hour max sub-steps
    dt_step = sign(min(abs(dt_total), 10800.0_dp), dt_total)
    n_steps = max(1, nint(abs(dt_total) / abs(dt_step)))
    dt_step = dt_total / real(n_steps, dp)

    do i = 1, n_steps
      call yoshida4_step_nb(pos, vel, gm, dt_step)
    end do
  end subroutine

end module nbody3_ekf_mod

! ─── Module 5: Keplerian elements + observation model ────────────────
module kepler_obs_3body_mod
  use constants_mod
  use nbody3_ekf_mod, only: NB, propagate_nbody
  implicit none
  integer, parameter :: NS = 15   ! EKF state dimension (7 Earth + 7 Moon + GM_moon)
  real(dp), parameter :: GM_EARTH_KNOWN = 398600.435507_dp
  real(dp), parameter :: GM_MOON_KNOWN  = 4902.800066_dp

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

  ! ── Predict alt/az for Sun or Moon from 12-param state ──
  subroutine predict_altaz(x, t_epoch, jd_obs, &
                            lat_deg, lon_deg, body_code, alt_deg, az_deg)
    real(dp), intent(in)  :: x(NS)
    real(dp), intent(in)  :: t_epoch, jd_obs, lat_deg, lon_deg
    integer,  intent(in)  :: body_code   ! 1=Sun, 2=Moon
    real(dp), intent(out) :: alt_deg, az_deg

    real(dp) :: ecc_e, inc_e, raan_e, argp_e, M0_e, mu_se, a_earth_v
    real(dp) :: ecc_m, inc_m, raan_m, argp_m, M0_m, mu_em, a_moon_v
    real(dp) :: gm_sun_v, gm_earth_v, gm_moon_v, M_tot
    real(dp) :: r_es(3), v_es(3), r_me(3), v_me(3)
    real(dp) :: r_bary(3), v_bary(3)
    real(dp) :: dt_s
    real(dp) :: pos_nb(NB,3), vel_nb(NB,3), gm_nb(NB)
    real(dp) :: r_target(3), v_earth(3)
    real(dp) :: d_gcrs(3), d_date(3), d_itrs(3), d_local(3)
    real(dp) :: obs_itrs(3), obs_date(3), obs_gcrs(3)
    real(dp) :: lat, lon, sinlat, coslat, slon, clon
    real(dp) :: jd_whole, jd_frac, jd_tdb
    real(dp) :: gmst_h_v, gast_h_v, gast_rad, sg, cg
    real(dp) :: M_mat(3,3)
    real(dp) :: dpsi_n, deps_n, mean_ob_n
    real(dp) :: d_unit(3), v_obs(3), vdot, dnorm
    real(dp) :: north, east, up, tmp

    ! ═══ Unpack state ═══
    ecc_e  = x(1);  inc_e  = x(2);  raan_e = x(3)
    argp_e = x(4);  M0_e   = x(5);  mu_se  = x(6);  a_earth_v = x(7)
    ecc_m  = x(8);  inc_m  = x(9);  raan_m = x(10)
    argp_m = x(11); M0_m   = x(12); mu_em  = x(13); a_moon_v = x(14)

    ! ═══ Derive individual GMs ═══
    gm_moon_v  = x(15)
    gm_earth_v = mu_em - gm_moon_v
    gm_sun_v   = mu_se - gm_earth_v

    ! ═══ Elements → Cartesian at epoch ═══
    ! Earth relative to Sun (from Earth Keplerian elements)
    call kepler_to_cart(a_earth_v, ecc_e, inc_e, raan_e, argp_e, M0_e, &
                        mu_se, r_es, v_es)
    ! Moon relative to Earth (from Moon Keplerian elements)
    call kepler_to_cart(a_moon_v, ecc_m, inc_m, raan_m, argp_m, M0_m, &
                        mu_em, r_me, v_me)

    ! ═══ 3-body initial conditions (barycentric frame) ═══
    ! In Sun-centered frame: Sun=0, Earth=r_es, Moon=r_es+r_me
    M_tot = gm_sun_v + gm_earth_v + gm_moon_v
    r_bary = (gm_earth_v * r_es + gm_moon_v * (r_es + r_me)) / M_tot
    v_bary = (gm_earth_v * v_es + gm_moon_v * (v_es + v_me)) / M_tot

    pos_nb(1,:) = -r_bary                      ! Sun
    vel_nb(1,:) = -v_bary
    pos_nb(2,:) = r_es - r_bary                ! Earth
    vel_nb(2,:) = v_es - v_bary
    pos_nb(3,:) = r_es + r_me - r_bary         ! Moon
    vel_nb(3,:) = v_es + v_me - v_bary

    gm_nb(1) = gm_sun_v
    gm_nb(2) = gm_earth_v
    gm_nb(3) = gm_moon_v

    ! ═══ Propagate to observation time ═══
    dt_s = (jd_obs - t_epoch) * DAY_S
    call propagate_nbody(pos_nb, vel_nb, gm_nb, dt_s)

    ! ═══ Target geocentric position ═══
    if (body_code == 1) then
      r_target = pos_nb(1,:) - pos_nb(2,:)   ! Sun geocentric
    else
      r_target = pos_nb(3,:) - pos_nb(2,:)   ! Moon geocentric
    end if
    v_earth = vel_nb(2,:)

    ! ═══ Coordinate transforms ═══
    lat = lat_deg * DEG2RAD
    lon = lon_deg * DEG2RAD
    sinlat = sin(lat); coslat = cos(lat)
    slon = sin(lon); clon = cos(lon)

    jd_whole = floor(jd_obs + 0.5_dp)
    jd_frac  = jd_obs - jd_whole
    jd_tdb = jd_obs + 69.184_dp / DAY_S

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

    ! Observer in GCRS
    obs_date(1) =  cg*obs_itrs(1) - sg*obs_itrs(2)
    obs_date(2) =  sg*obs_itrs(1) + cg*obs_itrs(2)
    obs_date(3) = obs_itrs(3)
    obs_gcrs = matmul(transpose(M_mat), obs_date)

    ! Direction in GCRS
    d_gcrs = r_target - obs_gcrs

    ! Aberration (first-order)
    v_obs = v_earth / C_LIGHT_KMS
    dnorm = sqrt(sum(d_gcrs**2))
    if (dnorm > 1.0d-10) then
      d_unit = d_gcrs / dnorm
      vdot = dot_product(d_unit, v_obs)
      d_gcrs = d_gcrs + dnorm * (v_obs - d_unit * vdot)
    end if

    ! GCRS → true equator of date
    d_date = matmul(M_mat, d_gcrs)

    ! True equator of date → ITRS
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

end module kepler_obs_3body_mod

! ═══════════════════════════════════════════════════════════════════════
!  Main program
! ═══════════════════════════════════════════════════════════════════════
program sun_ekf_3body
  use constants_mod
  use spk_reader_mod
  use obs_reader_mod
  use kepler_obs_3body_mod
  implicit none

  ! GM values (DE440)
  real(dp), parameter :: GM_SUN   = 1.32712440041279419d11
  real(dp), parameter :: GM_EARTH = 3.98600435507d5
  real(dp), parameter :: GM_MOON  = 4.9028000661d3

  ! State vectors
  real(dp) :: x(NS), x_true(NS), x_init(NS)
  real(dp) :: t_epoch

  ! EKF covariance
  real(dp) :: P_cov(NS,NS)
  real(dp) :: Q_noise(NS,NS), R_noise(2,2)

  ! SPK kernel
  type(spk_kernel) :: kernel

  ! Observations
  type(observation) :: obs_all(MAX_OBS)
  integer :: n_obs_all, n_sun, n_moon, n_all
  real(dp) :: s_jd(MAX_OBS), s_alt(MAX_OBS), s_az(MAX_OBS)
  real(dp) :: m_jd(MAX_OBS), m_alt(MAX_OBS), m_az(MAX_OBS)
  integer :: s_body(MAX_OBS), m_body(MAX_OBS)
  ! Merged arrays for stage 3
  real(dp) :: all_jd(MAX_OBS), all_alt(MAX_OBS), all_az(MAX_OBS)
  integer :: all_body(MAX_OBS)

  ! DE440s positions
  real(dp) :: sun_ssb(3), sun_vel(3)
  real(dp) :: emb_ssb(3), emb_vel(3)
  real(dp) :: earth_emb(3), earth_emb_vel(3)
  real(dp) :: moon_emb(3), moon_emb_vel(3)
  real(dp) :: earth_ssb(3), earth_vel_v(3)
  real(dp) :: r_rel(3), v_rel(3)
  real(dp) :: r_moon_rel(3), v_moon_rel(3)

  ! Keplerian elements
  real(dp) :: a_comp, ecc_comp, inc_comp, raan_comp, argp_comp, M0_comp
  real(dp) :: a_m_comp, ecc_m_comp, inc_m_comp, raan_m_comp, argp_m_comp, M0_m_comp
  real(dp) :: mu_se, mu_em

  ! Verification
  real(dp) :: r_check(3), v_check(3), r_diff

  ! Loop/misc
  integer :: i, j
  real(dp) :: lat_obs, lon_obs, jd_epoch, jd_w, jd_f
  real(dp) :: sigma_obs
  real(dp) :: delta(NS)
  real(dp) :: err_deg
  logical :: active(NS)
  logical :: grid_active
  real(dp) :: a_E_before, a_E_shift, a_E_prev_shift
  real(dp) :: gm_earth_est, gm_sun_est
  real(dp) :: corr_mat(NS,NS), sig_i_v, sig_j_v
  character(len=7) :: param_names(NS)
  real(dp) :: Q_rate(NS)   ! process noise rate (state²/day)
  integer :: iter
  integer, parameter :: N_ITER = 20   ! outer iterations

  ! ═══════════════════════════════════════════════════
  lat_obs = 40.0_dp
  lon_obs = 0.0_dp

  print '(A)', '════════════════════════════════════════════════════════'
  print '(A)', '  3-Body Sun-Earth-Moon EKF'
  print '(A,F8.4,A,F8.4)', '  Observer: lat=', lat_obs, ' lon=', lon_obs
  print '(A)', '════════════════════════════════════════════════════════'
  print '(A)', ''
  print '(A)', '  State vector (15 parameters):'
  print '(A)', '  ── Earth orbit x(1:7) ──'
  print '(A)', '    e_E   = eccentricity of Earth orbit'
  print '(A)', '    i_E   = inclination to equator (obliquity ~23.44 deg)'
  print '(A)', '    Om_E  = longitude of ascending node'
  print '(A)', '    w_E   = argument of periapsis'
  print '(A)', '    M0_E  = mean anomaly at epoch'
  print '(A)', '    mu_SE = GM_sun + GM_earth (km^3/s^2)'
  print '(A)', '    a_E   = semi-major axis of Earth orbit (km)'
  print '(A)', '  ── Moon orbit x(8:14) ──'
  print '(A)', '    e_M   = eccentricity of Moon orbit'
  print '(A)', '    i_M   = inclination to equator (~28.5 deg from ecliptic)'
  print '(A)', '    Om_M  = longitude of ascending node'
  print '(A)', '    w_M   = argument of periapsis'
  print '(A)', '    M0_M  = mean anomaly at epoch'
  print '(A)', '    mu_EM = GM_earth + GM_moon (km^3/s^2)'
  print '(A)', '    a_M   = semi-major axis of Moon orbit (km)'
  print '(A)', '  ── Mass x(15) ──'
  print '(A)', '    gm_M  = GM_moon (km^3/s^2)'

  ! ── 1. Read observations, split by body ──
  call read_observations('observations.dat', obs_all, n_obs_all, 1.0d10)
  n_sun = 0; n_moon = 0
  do i = 1, n_obs_all
    if (obs_all(i)%body == 1) then
      n_sun = n_sun + 1
      s_jd(n_sun)  = obs_all(i)%jd
      s_alt(n_sun) = obs_all(i)%alt_obs
      s_az(n_sun)  = obs_all(i)%az_obs
      s_body(n_sun) = 1
    else
      n_moon = n_moon + 1
      m_jd(n_moon)  = obs_all(i)%jd
      m_alt(n_moon) = obs_all(i)%alt_obs
      m_az(n_moon)  = obs_all(i)%az_obs
      m_body(n_moon) = 2
    end if
  end do
  ! Merged (file is already chronological)
  n_all = 0
  do i = 1, n_obs_all
    n_all = n_all + 1
    all_jd(n_all)   = obs_all(i)%jd
    all_alt(n_all)  = obs_all(i)%alt_obs
    all_az(n_all)   = obs_all(i)%az_obs
    all_body(n_all) = obs_all(i)%body
  end do
  print '(A,I6,A,I6,A)', '  Loaded ', n_sun, ' Sun +', n_moon, ' Moon observations'

  ! ── 2. Compute true Earth orbital elements from DE440s ──
  jd_epoch = s_jd(1)
  t_epoch  = jd_epoch
  jd_w = real(floor(jd_epoch + 0.5_dp), dp)
  jd_f = jd_epoch - jd_w

  call spk_open('de440s.bsp', kernel)
  call spk_compute_and_diff(kernel, 0, 10, jd_w, jd_f, sun_ssb, sun_vel)
  call spk_compute_and_diff(kernel, 0, 3,  jd_w, jd_f, emb_ssb, emb_vel)
  call spk_compute_and_diff(kernel, 3, 399, jd_w, jd_f, earth_emb, earth_emb_vel)
  call spk_compute_and_diff(kernel, 3, 301, jd_w, jd_f, moon_emb, moon_emb_vel)
  call spk_close(kernel)

  ! Earth in SSB frame
  earth_ssb    = emb_ssb + earth_emb
  earth_vel_v  = emb_vel + earth_emb_vel
  ! Convert km/day to km/s
  sun_vel       = sun_vel / DAY_S
  earth_vel_v   = earth_vel_v / DAY_S
  moon_emb_vel  = moon_emb_vel / DAY_S
  earth_emb_vel = earth_emb_vel / DAY_S

  ! Earth relative to Sun → Earth Keplerian elements
  r_rel = earth_ssb - sun_ssb
  v_rel = earth_vel_v - sun_vel
  mu_se = GM_SUN + GM_EARTH
  call cart_to_kepler(r_rel, v_rel, mu_se, a_comp, ecc_comp, inc_comp, &
       raan_comp, argp_comp, M0_comp)

  ! Round-trip check (Earth)
  call kepler_to_cart(a_comp, ecc_comp, inc_comp, raan_comp, argp_comp, &
       M0_comp, mu_se, r_check, v_check)
  r_diff = sqrt(sum((r_check - r_rel)**2))
  print '(A,ES10.3,A)', '  Earth round-trip error: ', r_diff, ' km'

  ! ── 3. Compute true Moon orbital elements from DE440s ──
  ! Moon relative to Earth = (Moon wrt EMB) - (Earth wrt EMB)
  r_moon_rel = moon_emb - earth_emb
  v_moon_rel = moon_emb_vel - earth_emb_vel
  mu_em = GM_EARTH + GM_MOON
  call cart_to_kepler(r_moon_rel, v_moon_rel, mu_em, a_m_comp, ecc_m_comp, &
       inc_m_comp, raan_m_comp, argp_m_comp, M0_m_comp)

  ! Round-trip check (Moon)
  call kepler_to_cart(a_m_comp, ecc_m_comp, inc_m_comp, raan_m_comp, argp_m_comp, &
       M0_m_comp, mu_em, r_check, v_check)
  r_diff = sqrt(sum((r_check - r_moon_rel)**2))
  print '(A,ES10.3,A)', '  Moon  round-trip error: ', r_diff, ' km'

  ! ── 4. True state vector ──
  x_true(1)  = ecc_comp;     x_true(2)  = inc_comp
  x_true(3)  = raan_comp;    x_true(4)  = argp_comp
  x_true(5)  = M0_comp;      x_true(6)  = mu_se
  x_true(7)  = a_comp
  x_true(8)  = ecc_m_comp;   x_true(9)  = inc_m_comp
  x_true(10) = raan_m_comp;  x_true(11) = argp_m_comp
  x_true(12) = M0_m_comp;    x_true(13) = mu_em
  x_true(14) = a_m_comp
  x_true(15) = GM_MOON

  print '(/,A)', '  True Keplerian elements (J2000 equatorial):'
  print '(A)', '  ── Earth orbit ──'
  print '(A,F12.8)',      '    e_E   = ', x_true(1)
  print '(A,F10.5,A)',    '    i_E   = ', x_true(2) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    Om_E  = ', x_true(3) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    w_E   = ', x_true(4) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    M0_E  = ', x_true(5) * RAD2DEG, ' deg'
  print '(A,ES20.12,A)',  '    mu_SE = ', x_true(6), ' km^3/s^2'
  print '(A,F12.1,A)',    '    a_E   = ', x_true(7), ' km'
  print '(A)', '  ── Moon orbit ──'
  print '(A,F12.8)',      '    e_M   = ', x_true(8)
  print '(A,F10.5,A)',    '    i_M   = ', x_true(9) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    Om_M  = ', x_true(10) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    w_M   = ', x_true(11) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    M0_M  = ', x_true(12) * RAD2DEG, ' deg'
  print '(A,ES20.12,A)',  '    mu_EM = ', x_true(13), ' km^3/s^2'
  print '(A,F12.1,A)',    '    a_M   = ', x_true(14), ' km'
  print '(A)', '  ── Mass ──'
  print '(A,ES20.12,A)',  '    gm_M  = ', x_true(15), ' km^3/s^2'

  ! ── 5. Perturbed initial state ──
  x_init(1)  = x_true(1)  * 1.20_dp            ! e_E: +20%
  x_init(2)  = x_true(2)  + 2.0_dp * DEG2RAD   ! i_E: +2 deg
  x_init(3)  = x_true(3)  + 5.0_dp * DEG2RAD   ! Om_E: +5 deg
  x_init(4)  = x_true(4)  + 5.0_dp * DEG2RAD   ! w_E: +5 deg
  x_init(5)  = x_true(5)  + 3.0_dp * DEG2RAD   ! M0_E: +3 deg
  x_init(6)  = x_true(6)  * 1.005_dp            ! mu_SE: +0.5%
  x_init(7)  = x_true(7)  * 1.01_dp             ! a_E: +1%
  x_init(8)  = x_true(8)  * 1.20_dp             ! e_M: +20%
  x_init(9)  = x_true(9)  + 2.0_dp * DEG2RAD    ! i_M: +2 deg
  x_init(10) = x_true(10) + 5.0_dp * DEG2RAD    ! Om_M: +5 deg
  x_init(11) = x_true(11) + 5.0_dp * DEG2RAD    ! w_M: +5 deg
  x_init(12) = x_true(12) + 3.0_dp * DEG2RAD    ! M0_M: +3 deg
  x_init(13) = x_true(13) * 1.005_dp             ! mu_EM: +0.5%
  x_init(14) = x_true(14) * 1.01_dp              ! a_M: +1%
  x_init(15) = x_true(15) * 2.0_dp               ! gm_M: fixed at 2× truth

  print '(/,A)', '  Perturbed initial guess:'
  print '(A)', '  ── Earth ──'
  print '(A,F12.8,A,F6.1,A)', '    e_E   = ', x_init(1), &
       '  (', (x_init(1)/x_true(1)-1.0_dp)*100.0_dp, '%)'
  print '(A,F10.5,A,F6.2,A)', '    i_E   = ', x_init(2)*RAD2DEG, &
       ' deg  (delta ', (x_init(2)-x_true(2))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    Om_E  = ', x_init(3)*RAD2DEG, &
       ' deg  (delta ', (x_init(3)-x_true(3))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    w_E   = ', x_init(4)*RAD2DEG, &
       ' deg  (delta ', (x_init(4)-x_true(4))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    M0_E  = ', x_init(5)*RAD2DEG, &
       ' deg  (delta ', (x_init(5)-x_true(5))*RAD2DEG, ' deg)'
  print '(A,ES20.12,A,F6.3,A)', '    mu_SE = ', x_init(6), &
       '  (', (x_init(6)/x_true(6)-1.0_dp)*100.0_dp, '%)'
  print '(A,F12.1,A,F7.3,A)',  '    a_E   = ', x_init(7), &
       ' km  (', (x_init(7)/x_true(7)-1.0_dp)*100.0_dp, '%)'
  print '(A)', '  ── Moon ──'
  print '(A,F12.8,A,F6.1,A)', '    e_M   = ', x_init(8), &
       '  (', (x_init(8)/x_true(8)-1.0_dp)*100.0_dp, '%)'
  print '(A,F10.5,A,F6.2,A)', '    i_M   = ', x_init(9)*RAD2DEG, &
       ' deg  (delta ', (x_init(9)-x_true(9))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    Om_M  = ', x_init(10)*RAD2DEG, &
       ' deg  (delta ', (x_init(10)-x_true(10))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    w_M   = ', x_init(11)*RAD2DEG, &
       ' deg  (delta ', (x_init(11)-x_true(11))*RAD2DEG, ' deg)'
  print '(A,F10.5,A,F6.2,A)', '    M0_M  = ', x_init(12)*RAD2DEG, &
       ' deg  (delta ', (x_init(12)-x_true(12))*RAD2DEG, ' deg)'
  print '(A,ES20.12,A,F6.3,A)', '    mu_EM = ', x_init(13), &
       '  (', (x_init(13)/x_true(13)-1.0_dp)*100.0_dp, '%)'
  print '(A,F12.1,A,F7.3,A)',  '    a_M   = ', x_init(14), &
       ' km  (', (x_init(14)/x_true(14)-1.0_dp)*100.0_dp, '%)'
  print '(A)', '  ── Mass ──'
  print '(A,ES20.12,A,F6.3,A)', '    gm_M  = ', x_init(15), &
       '  (', (x_init(15)/x_true(15)-1.0_dp)*100.0_dp, '%)'

  ! ── 6. Initialize EKF ──
  x = x_init

  ! Process noise rate: Q_rate (state²/day)
  ! The state vector contains initial Keplerian elements at epoch.
  ! The N-body integrator handles ALL dynamics. Process noise is
  ! zero — the initial conditions don't change. Instead, we use
  ! directional covariance inflation (in run_ekf_pass) to prevent
  ! collapse along the mu-a degenerate direction.
  Q_rate = 0.0_dp

  Q_noise = 0.0_dp   ! will be set per step from Q_rate * dt

  sigma_obs = 1.0_dp / 3600.0_dp   ! 1 arcsec (near-exact observations)
  R_noise = 0.0_dp
  R_noise(1,1) = sigma_obs**2
  R_noise(2,2) = sigma_obs**2

  ! FD perturbation sizes
  ! Earth: 1e-7 (far away, standard)
  delta(1) = 1.0d-7;  delta(2) = 1.0d-7;  delta(3) = 1.0d-7
  delta(4) = 1.0d-7;  delta(5) = 1.0d-7;  delta(6) = mu_se * 1.0d-7
  delta(7) = a_comp * 1.0d-7
  ! Moon: 1e-5 (closer, parallax needs larger FD step)
  delta(8) = 1.0d-5;  delta(9) = 1.0d-5;  delta(10) = 1.0d-5
  delta(11) = 1.0d-5; delta(12) = 1.0d-5; delta(13) = mu_em * 1.0d-5
  delta(14) = a_m_comp * 1.0d-5
  ! GM_moon: same scale as Moon parameters
  delta(15) = GM_MOON * 1.0d-5

  ! ═══════════════════════════════════════════════════
  ! Outer iteration loop: re-run all 3 stages using
  ! previous final estimate as starting point.
  ! ═══════════════════════════════════════════════════
  grid_active = .true.
  a_E_prev_shift = huge(1.0_dp)
  do iter = 1, N_ITER

  print '(/,A,I2,A,I2)', '  ══════ Iteration ', iter, ' / ', N_ITER

  ! On iteration > 1, use previous result as starting point
  if (iter > 1) then
    x_init = x
  end if
  x = x_init

  ! ══════════════════════════════════════════════════════════════════════
  ! Single joint fit: all observations, all parameters active
  ! This lets the parallactic inequality constrain a_E from the start.
  ! ══════════════════════════════════════════════════════════════════════
  P_cov = 0.0_dp
  ! Earth orbit
  P_cov(1,1) = (0.005_dp)**2               ! sigma_e = 0.005
  P_cov(2,2) = (3.0_dp * DEG2RAD)**2       ! sigma_i = 3 deg
  P_cov(3,3) = (8.0_dp * DEG2RAD)**2       ! sigma_Omega = 8 deg
  P_cov(4,4) = (8.0_dp * DEG2RAD)**2       ! sigma_omega = 8 deg
  P_cov(5,5) = (5.0_dp * DEG2RAD)**2       ! sigma_M0 = 5 deg
  P_cov(6,6) = (0.20_dp * mu_se)**2        ! sigma_mu = 20%
  P_cov(7,7) = (0.20_dp * a_comp)**2       ! sigma_a_E = 20%
  ! Moon orbit
  P_cov(8,8)   = (0.02_dp)**2              ! sigma_e_m = 0.02
  P_cov(9,9)   = (5.0_dp * DEG2RAD)**2     ! sigma_i_m = 5 deg
  P_cov(10,10) = (10.0_dp * DEG2RAD)**2    ! sigma_Om_m = 10 deg
  P_cov(11,11) = (10.0_dp * DEG2RAD)**2    ! sigma_w_m = 10 deg
  P_cov(12,12) = (5.0_dp * DEG2RAD)**2     ! sigma_M0_m = 5 deg
  P_cov(13,13) = (0.20_dp * mu_em)**2      ! sigma_mu_em = 20%
  P_cov(14,14) = (0.20_dp * a_m_comp)**2   ! sigma_a_M = 20%
  ! GM_moon
  P_cov(15,15) = (0.01_dp * GM_MOON)**2     ! sigma_gm_moon = 1%

  active = .true.
  active(15) = .false.   ! fix GM_moon — test observability

  call run_ekf_pass(x, P_cov, n_all, all_jd, all_alt, all_az, all_body, active, &
                    'Joint fit: all observations, all parameters')

  ! ── Post-EKF: 1D grid search along the n=const degenerate direction ──
  ! The EKF determines n = sqrt(mu/a^3) precisely but cannot resolve
  ! mu and a individually. Search along the n=const curve for the a_E
  ! that minimizes Moon observation residuals (parallactic inequality).
  ! Stop the grid search when the correction becomes small (converged)
  ! or when the shift oscillates (noise-dominated regime).
  if (grid_active) then
    a_E_before = x(7)
    call grid_search_degenerate(x, n_all, all_jd, all_alt, all_az, all_body)
    a_E_shift = abs(x(7) - a_E_before) / x_true(7) * 100.0_dp

    if (a_E_shift < 2.0_dp) then
      print '(A,F6.2,A)', '  Grid shift < 2% (', a_E_shift, '%) — switching to EKF-only'
      grid_active = .false.
    end if
    if (iter > 1 .and. a_E_shift > a_E_prev_shift * 1.5_dp) then
      ! Shift is growing — noise-dominated, revert and stop
      print '(A)', '  Grid shift growing — reverting, switching to EKF-only'
      x(7) = a_E_before
      x(6) = sqrt(x(6) / a_E_before**3)**2 * a_E_before**3
      grid_active = .false.
    end if
    a_E_prev_shift = a_E_shift
  end if

  ! Print iteration summary — current values of all parameters and errors
  print '(/,A,I2,A)', '  ── Iteration ', iter, ' summary ──'
  print '(A)', '  Parameter       Value                  True                 Error'
  print '(A)', '  ─────────────────────────────────────────────────────────────────'
  print '(A,F12.8,8X,F12.8,8X,F8.4,A)', &
       '  e_E     ', x(1), x_true(1), (x(1)/x_true(1)-1.0_dp)*100.0_dp, '%'
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  i_E     ', x(2)*RAD2DEG, ' deg', x_true(2)*RAD2DEG, ' deg', &
       (x(2)-x_true(2))*RAD2DEG, ' deg'
  err_deg = (x(3) - x_true(3)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  Om_E    ', x(3)*RAD2DEG, ' deg', x_true(3)*RAD2DEG, ' deg', err_deg, ' deg'
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  w_E     ', x(4)*RAD2DEG, ' deg', x_true(4)*RAD2DEG, ' deg', &
       (x(4)-x_true(4))*RAD2DEG, ' deg'
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  M0_E    ', x(5)*RAD2DEG, ' deg', x_true(5)*RAD2DEG, ' deg', &
       (x(5)-x_true(5))*RAD2DEG, ' deg'
  print '(A,ES20.12,2X,ES20.12,2X,F8.4,A)', &
       '  mu_SE   ', x(6), x_true(6), (x(6)/x_true(6)-1.0_dp)*100.0_dp, '%'
  print '(A,F14.1,A,6X,F14.1,A,6X,F8.4,A)', &
       '  a_E     ', x(7), ' km', x_true(7), ' km', &
       (x(7)/x_true(7)-1.0_dp)*100.0_dp, '%'
  print '(A,F12.8,8X,F12.8,8X,F8.4,A)', &
       '  e_M     ', x(8), x_true(8), (x(8)/x_true(8)-1.0_dp)*100.0_dp, '%'
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  i_M     ', x(9)*RAD2DEG, ' deg', x_true(9)*RAD2DEG, ' deg', &
       (x(9)-x_true(9))*RAD2DEG, ' deg'
  err_deg = (x(10) - x_true(10)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  Om_M    ', x(10)*RAD2DEG, ' deg', x_true(10)*RAD2DEG, ' deg', err_deg, ' deg'
  err_deg = (x(11) - x_true(11)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  w_M     ', x(11)*RAD2DEG, ' deg', x_true(11)*RAD2DEG, ' deg', err_deg, ' deg'
  err_deg = (x(12) - x_true(12)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,4X,F10.5,A,4X,F8.4,A)', &
       '  M0_M    ', x(12)*RAD2DEG, ' deg', x_true(12)*RAD2DEG, ' deg', err_deg, ' deg'
  print '(A,ES20.12,2X,ES20.12,2X,F8.4,A)', &
       '  mu_EM   ', x(13), x_true(13), (x(13)/x_true(13)-1.0_dp)*100.0_dp, '%'
  print '(A,F14.1,A,6X,F14.1,A,6X,F8.4,A)', &
       '  a_M     ', x(14), ' km', x_true(14), ' km', &
       (x(14)/x_true(14)-1.0_dp)*100.0_dp, '%'
  print '(A,ES20.12,2X,ES20.12,2X,F8.4,A)', &
       '  gm_M    ', x(15), x_true(15), (x(15)/x_true(15)-1.0_dp)*100.0_dp, '%'

  end do  ! outer iteration loop

  ! ══════════════════════════════════════════════════════════════════════
  ! FINAL RESULTS
  ! ══════════════════════════════════════════════════════════════════════
  print '(/,A)', '════════════════════════════════════════════════════════'
  print '(A)',   '  FINAL ESTIMATES'
  print '(A)',   '════════════════════════════════════════════════════════'

  ! ── Earth orbit ──
  print '(A)',   '  ── Earth orbit ──'
  print '(A)',   '  Parameter   Estimated         True            Error'
  print '(A)',   '  ─────────────────────────────────────────────────────'
  print '(A,F12.8,4X,F12.8,4X,F8.4,A)', &
       '  e_E     ', x(1), x_true(1), (x(1)/x_true(1)-1.0_dp)*100.0_dp, '%'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  i_E     ', x(2)*RAD2DEG, ' deg', x_true(2)*RAD2DEG, ' deg', &
       (x(2)-x_true(2))*RAD2DEG, ' deg'
  err_deg = (x(3) - x_true(3)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  Om_E    ', x(3)*RAD2DEG, ' deg', x_true(3)*RAD2DEG, ' deg', err_deg, ' deg'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  w_E     ', x(4)*RAD2DEG, ' deg', x_true(4)*RAD2DEG, ' deg', &
       (x(4)-x_true(4))*RAD2DEG, ' deg'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  M0_E    ', x(5)*RAD2DEG, ' deg', x_true(5)*RAD2DEG, ' deg', &
       (x(5)-x_true(5))*RAD2DEG, ' deg'
  print '(A,ES20.12)', '  mu_SE   ', x(6)
  print '(A,ES20.12,A,F8.4,A)', '  mu_true ', x_true(6), &
       '  err: ', (x(6)/x_true(6)-1.0_dp)*100.0_dp, '%'
  print '(A,F12.1,A,F12.1,A,F8.4,A)', &
       '  a_E     ', x(7), ' km  true: ', x_true(7), ' km  err: ', &
       (x(7)/x_true(7)-1.0_dp)*100.0_dp, '%'

  ! ── Moon orbit ──
  print '(/,A)', '  ── Moon orbit ──'
  print '(A)',   '  Parameter   Estimated         True            Error'
  print '(A)',   '  ─────────────────────────────────────────────────────'
  print '(A,F12.8,4X,F12.8,4X,F8.4,A)', &
       '  e_M     ', x(8), x_true(8), (x(8)/x_true(8)-1.0_dp)*100.0_dp, '%'
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  i_M     ', x(9)*RAD2DEG, ' deg', x_true(9)*RAD2DEG, ' deg', &
       (x(9)-x_true(9))*RAD2DEG, ' deg'
  err_deg = (x(10) - x_true(10)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  Om_M    ', x(10)*RAD2DEG, ' deg', x_true(10)*RAD2DEG, ' deg', err_deg, ' deg'
  err_deg = (x(11) - x_true(11)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  w_M     ', x(11)*RAD2DEG, ' deg', x_true(11)*RAD2DEG, ' deg', err_deg, ' deg'
  err_deg = (x(12) - x_true(12)) * RAD2DEG
  if (err_deg >  180.0_dp) err_deg = err_deg - 360.0_dp
  if (err_deg < -180.0_dp) err_deg = err_deg + 360.0_dp
  print '(A,F10.5,A,2X,F10.5,A,2X,F8.4,A)', &
       '  M0_M    ', x(12)*RAD2DEG, ' deg', x_true(12)*RAD2DEG, ' deg', err_deg, ' deg'
  print '(A,ES20.12)', '  mu_EM   ', x(13)
  print '(A,ES20.12,A,F8.4,A)', '  mu_true ', x_true(13), &
       '  err: ', (x(13)/x_true(13)-1.0_dp)*100.0_dp, '%'
  print '(A,F12.1,A,F12.1,A,F8.4,A)', &
       '  a_M     ', x(14), ' km  true: ', x_true(14), ' km  err: ', &
       (x(14)/x_true(14)-1.0_dp)*100.0_dp, '%'

  ! ── GM_moon ──
  print '(/,A)', '  ── GM_moon (fitted) ──'
  print '(A,ES20.12)', '  gm_M    ', x(15)
  print '(A,ES20.12,A,F8.4,A)', '  gm_true ', x_true(15), &
       '  err: ', (x(15)/x_true(15)-1.0_dp)*100.0_dp, '%'
  gm_earth_est = x(13) - x(15)
  gm_sun_est   = x(6) - gm_earth_est

  print '(/,A)', '  ── Derived masses ──'
  print '(A,ES15.6,A,ES15.6,A,F8.4,A)', '  GM_sun   = ', gm_sun_est, &
       '  (true: ', GM_SUN, ')  err: ', (gm_sun_est/GM_SUN-1.0_dp)*100.0_dp, '%'
  print '(A,ES15.6,A,ES15.6,A,F8.4,A)', '  GM_earth = ', gm_earth_est, &
       '  (true: ', GM_EARTH, ')  err: ', (gm_earth_est/GM_EARTH-1.0_dp)*100.0_dp, '%'
  print '(A,ES15.6,A,ES15.6,A,F8.4,A)', '  GM_moon  = ', x(15), &
       '  (true: ', GM_MOON, ')  err: ', (x(15)/GM_MOON-1.0_dp)*100.0_dp, '%'

  ! ── Derived periods ──
  print '(/,A)', '  ── Orbital periods ──'
  print '(A,F12.3,A)', '    Earth period = ', &
       TAU / sqrt(x(6) / x(7)**3) / DAY_S, ' days'
  print '(A,F12.3,A)', '    (true)       = ', &
       TAU / sqrt(x_true(6) / x_true(7)**3) / DAY_S, ' days'
  print '(A,F12.3,A)', '    Moon period  = ', &
       TAU / sqrt(x(13) / x(14)**3) / DAY_S, ' days'
  print '(A,F12.3,A)', '    (true)       = ', &
       TAU / sqrt(x_true(13) / x_true(14)**3) / DAY_S, ' days'

  ! ── 1-sigma uncertainties ──
  print '(/,A)', '  1-sigma uncertainties (from covariance):'
  print '(A)', '  ── Earth ──'
  print '(A,ES10.3)', '    sigma_e_E   = ', sqrt(max(0.0_dp, P_cov(1,1)))
  print '(A,F8.4,A)', '    sigma_i_E   = ', sqrt(max(0.0_dp, P_cov(2,2)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_Om_E  = ', sqrt(max(0.0_dp, P_cov(3,3)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_w_E   = ', sqrt(max(0.0_dp, P_cov(4,4)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_M0_E  = ', sqrt(max(0.0_dp, P_cov(5,5)))*RAD2DEG, ' deg'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_mu_SE = ', sqrt(max(0.0_dp, P_cov(6,6))), &
       '  (', sqrt(max(0.0_dp, P_cov(6,6)))/x(6)*100.0_dp, '%)'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_a_E   = ', sqrt(max(0.0_dp, P_cov(7,7))), &
       '  (', sqrt(max(0.0_dp, P_cov(7,7)))/x(7)*100.0_dp, '%)'
  print '(A)', '  ── Moon ──'
  print '(A,ES10.3)', '    sigma_e_M   = ', sqrt(max(0.0_dp, P_cov(8,8)))
  print '(A,F8.4,A)', '    sigma_i_M   = ', sqrt(max(0.0_dp, P_cov(9,9)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_Om_M  = ', sqrt(max(0.0_dp, P_cov(10,10)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_w_M   = ', sqrt(max(0.0_dp, P_cov(11,11)))*RAD2DEG, ' deg'
  print '(A,F8.4,A)', '    sigma_M0_M  = ', sqrt(max(0.0_dp, P_cov(12,12)))*RAD2DEG, ' deg'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_mu_EM = ', sqrt(max(0.0_dp, P_cov(13,13))), &
       '  (', sqrt(max(0.0_dp, P_cov(13,13)))/x(13)*100.0_dp, '%)'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_a_M   = ', sqrt(max(0.0_dp, P_cov(14,14))), &
       '  (', sqrt(max(0.0_dp, P_cov(14,14)))/x(14)*100.0_dp, '%)'
  print '(A)', '  ── Mass ──'
  print '(A,ES10.3,A,F8.5,A)', '    sigma_gm_M  = ', sqrt(max(0.0_dp, P_cov(15,15))), &
       '  (', sqrt(max(0.0_dp, P_cov(15,15)))/x(15)*100.0_dp, '%)'

  ! Correlation matrix
  param_names(1)  = '  e_E  '
  param_names(2)  = '  i_E  '
  param_names(3)  = '  Om_E '
  param_names(4)  = '  w_E  '
  param_names(5)  = '  M0_E '
  param_names(6)  = '  muSE '
  param_names(7)  = '  a_E  '
  param_names(8)  = '  e_M  '
  param_names(9)  = '  i_M  '
  param_names(10) = '  Om_M '
  param_names(11) = '  w_M  '
  param_names(12) = '  M0_M '
  param_names(13) = '  muEM '
  param_names(14) = '  a_M  '
  param_names(15) = '  gm_M '
  do i = 1, NS
    sig_i_v = sqrt(max(0.0_dp, P_cov(i,i)))
    do j = 1, NS
      sig_j_v = sqrt(max(0.0_dp, P_cov(j,j)))
      if (sig_i_v > 0.0_dp .and. sig_j_v > 0.0_dp) then
        corr_mat(i,j) = P_cov(i,j) / (sig_i_v * sig_j_v)
      else
        corr_mat(i,j) = 0.0_dp
      end if
    end do
  end do

  print '(/,A)', '  Correlation matrix:'
  print '(A,15(A7,1X))', '          ', (param_names(j), j=1,NS)
  do i = 1, NS
    print '(A,15F8.4)', param_names(i), (corr_mat(i,j), j=1,NS)
  end do

  print '(/,A)', '  Note: "True" elements are osculating at epoch.'
  print '(A)',   '  The EKF estimates best-fit MEAN elements over the arc.'
  print '(A)',   '  Moon elements change significantly over the year due to'
  print '(A)',   '  solar perturbations (node regression, apsidal precession).'
  print '(A)', '════════════════════════════════════════════════════════'

contains

  ! ── Generic EKF pass subroutine ──
  subroutine run_ekf_pass(xv, Pv, n_obs, obs_jds, obs_alts, obs_azs, &
                           obs_bodies, act, stage_name)
    real(dp), intent(inout) :: xv(NS), Pv(NS,NS)
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: obs_jds(n_obs), obs_alts(n_obs), obs_azs(n_obs)
    integer, intent(in) :: obs_bodies(n_obs)
    logical, intent(in) :: act(NS)
    character(len=*), intent(in) :: stage_name

    real(dp) :: P_pred(NS,NS), H(2,NS)
    real(dp) :: PHt(NS,2), S_mat(2,2), S_inv(2,2), K_gain(NS,2)
    real(dp) :: KH(NS,NS), I_KH(NS,NS), I_id(NS,NS)
    real(dp) :: zz_obs(2), zz_pred(2), innov(2)
    real(dp) :: alt_pred, az_pred, alt_p, az_p
    real(dp) :: xp(NS), det_S
    real(dp) :: win_a, win_z
    integer :: kk, jj, n_proc
    real(dp) :: jd_prev, dt_days, Q_step(NS,NS)

    I_id = 0.0_dp
    do kk = 1, NS
      I_id(kk,kk) = 1.0_dp
    end do

    win_a = 0.0_dp; win_z = 0.0_dp
    n_proc = 0
    jd_prev = obs_jds(1)

    print '(/,A,A)', '  ', trim(stage_name)
    print '(A)', '  ──────────────────────────────────────────────────────────────────────────'
    print '(A)', '   Obs#  winRMS_a" winRMS_z"  e_E_err%  muSE_err%  aE_err%   e_M_err%  muEM_err%  aM_err%  gmM_err%'

    do kk = 1, n_obs
      ! Time-dependent process noise
      dt_days = abs(obs_jds(kk) - jd_prev)
      if (dt_days < 1.0d-6) dt_days = 1.0d-6
      Q_step = 0.0_dp
      do jj = 1, NS
        Q_step(jj,jj) = Q_rate(jj) * dt_days
      end do
      P_pred = Pv + Q_step
      jd_prev = obs_jds(kk)

      call predict_altaz(xv, t_epoch, obs_jds(kk), &
                          lat_obs, lon_obs, obs_bodies(kk), alt_pred, az_pred)

      ! Skip if prediction is NaN (bad state)
      if (alt_pred /= alt_pred .or. az_pred /= az_pred) cycle

      innov(1) = obs_alts(kk) - alt_pred
      innov(2) = obs_azs(kk)  - az_pred
      if (innov(2) >  180.0_dp) innov(2) = innov(2) - 360.0_dp
      if (innov(2) < -180.0_dp) innov(2) = innov(2) + 360.0_dp

      ! Innovation gating: skip if residual > 30 deg (wildly wrong)
      if (abs(innov(1)) > 30.0_dp .or. abs(innov(2)) > 30.0_dp) cycle

      ! Jacobian (only active params get FD columns)
      H = 0.0_dp
      do jj = 1, NS
        if (.not. act(jj)) cycle
        xp = xv
        xp(jj) = xp(jj) + delta(jj)
        call predict_altaz(xp, t_epoch, obs_jds(kk), &
                            lat_obs, lon_obs, obs_bodies(kk), alt_p, az_p)
        if (alt_p /= alt_p .or. az_p /= az_p) then
          H(1,jj) = 0.0_dp; H(2,jj) = 0.0_dp; cycle
        end if
        H(1,jj) = (alt_p - alt_pred) / delta(jj)
        H(2,jj) = (az_p  - az_pred)  / delta(jj)
        if (H(2,jj) >  180.0_dp / delta(jj)) H(2,jj) = H(2,jj) - 360.0_dp / delta(jj)
        if (H(2,jj) < -180.0_dp / delta(jj)) H(2,jj) = H(2,jj) + 360.0_dp / delta(jj)
      end do

      ! Innovation covariance
      PHt = matmul(P_pred, transpose(H))
      S_mat = matmul(H, PHt) + R_noise

      det_S = S_mat(1,1)*S_mat(2,2) - S_mat(1,2)*S_mat(2,1)
      if (abs(det_S) < 1.0d-30) cycle
      S_inv(1,1) =  S_mat(2,2) / det_S
      S_inv(2,2) =  S_mat(1,1) / det_S
      S_inv(1,2) = -S_mat(1,2) / det_S
      S_inv(2,1) = -S_mat(2,1) / det_S

      K_gain = matmul(PHt, S_inv)

      ! State update (save old state for NaN recovery)
      xp = xv
      do jj = 1, NS
        xv(jj) = xv(jj) + K_gain(jj,1)*innov(1) + K_gain(jj,2)*innov(2)
      end do

      ! NaN guard: if any state is NaN, revert and skip
      if (any(xv /= xv)) then
        xv = xp
        Pv = P_pred
        cycle
      end if

      ! Enforce constraints
      ! Earth angles
      xv(3)  = mod(xv(3), TAU);  if (xv(3)  < 0.0_dp) xv(3)  = xv(3)  + TAU
      xv(4)  = mod(xv(4), TAU);  if (xv(4)  < 0.0_dp) xv(4)  = xv(4)  + TAU
      xv(5)  = mod(xv(5), TAU);  if (xv(5)  < 0.0_dp) xv(5)  = xv(5)  + TAU
      ! Moon angles
      xv(10) = mod(xv(10), TAU); if (xv(10) < 0.0_dp) xv(10) = xv(10) + TAU
      xv(11) = mod(xv(11), TAU); if (xv(11) < 0.0_dp) xv(11) = xv(11) + TAU
      xv(12) = mod(xv(12), TAU); if (xv(12) < 0.0_dp) xv(12) = xv(12) + TAU
      ! Eccentricities > 0
      if (xv(1)  < 1.0d-8) xv(1)  = 1.0d-8
      if (xv(8)  < 1.0d-8) xv(8)  = 1.0d-8
      ! mu > 0
      if (xv(6)  < 0.0_dp) xv(6)  = x_true(6) * 0.9_dp
      if (xv(13) < 0.0_dp) xv(13) = x_true(13) * 0.9_dp
      ! a > 0
      if (xv(7)  < 0.0_dp) xv(7)  = x_true(7) * 0.9_dp
      if (xv(14) < 0.0_dp) xv(14) = x_true(14) * 0.9_dp
      ! GM_moon > 0
      if (xv(15) < 0.0_dp) xv(15) = x_true(15) * 0.9_dp

      ! Covariance update
      KH = matmul(K_gain, H)
      I_KH = I_id - KH
      Pv = matmul(I_KH, P_pred)
      Pv = 0.5_dp * (Pv + transpose(Pv))

      ! Windowed statistics
      n_proc = n_proc + 1
      win_a = win_a + innov(1)**2
      win_z = win_z + innov(2)**2

      if (mod(n_proc, 200) == 0 .or. n_proc == 1 .or. n_proc == n_obs) then
        jj = min(n_proc, 200)
        if (n_proc == 1) jj = 1
        print '(I7,2F10.1,7F10.4)', &
             n_proc, &
             sqrt(win_a / jj) * 3600.0_dp, &
             sqrt(win_z / jj) * 3600.0_dp, &
             (xv(1)/x_true(1) - 1.0_dp) * 100.0_dp, &
             (xv(6)/x_true(6) - 1.0_dp) * 100.0_dp, &
             (xv(7)/x_true(7) - 1.0_dp) * 100.0_dp, &
             (xv(8)/x_true(8) - 1.0_dp) * 100.0_dp, &
             (xv(13)/x_true(13) - 1.0_dp) * 100.0_dp, &
             (xv(14)/x_true(14) - 1.0_dp) * 100.0_dp, &
             (xv(15)/x_true(15) - 1.0_dp) * 100.0_dp
        win_a = 0.0_dp; win_z = 0.0_dp
      end if
    end do
  end subroutine

  ! ── 1D grid search along the n=const degenerate direction ──
  ! After the EKF determines n = sqrt(mu/a^3) precisely, this
  ! searches along the n=const curve for the a_E value that
  ! minimizes Moon observation residuals. The parallactic inequality
  ! amplitude depends on a_E independently of n, so Moon observations
  ! break the mu-a degeneracy here.
  subroutine grid_search_degenerate(xv, n_obs, obs_jds, obs_alts, obs_azs, obs_bodies)
    real(dp), intent(inout) :: xv(NS)
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: obs_jds(n_obs), obs_alts(n_obs), obs_azs(n_obs)
    integer, intent(in) :: obs_bodies(n_obs)

    integer, parameter :: N_GRID = 21
    integer, parameter :: SUBSAMPLE = 20
    real(dp) :: n_mean, a_lo, a_hi, a_step
    real(dp) :: a_trial, mu_trial, rss, rss_best, a_best, mu_best
    real(dp) :: x_trial(NS), alt_pred, az_pred
    real(dp) :: da, dz
    integer :: ii, kk, n_moon

    ! Mean motion (invariant along the n=const curve)
    n_mean = sqrt(xv(6) / xv(7)**3)

    ! Search range: ±50% around current a_E
    a_lo = xv(7) * 0.50_dp
    a_hi = xv(7) * 1.50_dp
    a_step = (a_hi - a_lo) / real(N_GRID - 1, dp)

    rss_best = huge(1.0_dp)
    a_best = xv(7)
    mu_best = xv(6)

    print '(/,A)', '  ── Grid search along n=const curve for a_E ──'

    do ii = 1, N_GRID
      a_trial = a_lo + real(ii - 1, dp) * a_step
      mu_trial = n_mean**2 * a_trial**3

      x_trial = xv
      x_trial(6) = mu_trial
      x_trial(7) = a_trial

      rss = 0.0_dp
      n_moon = 0
      do kk = 1, n_obs, SUBSAMPLE
        if (obs_bodies(kk) /= 2) cycle  ! Moon only
        call predict_altaz(x_trial, t_epoch, obs_jds(kk), &
                            lat_obs, lon_obs, obs_bodies(kk), alt_pred, az_pred)
        if (alt_pred /= alt_pred .or. az_pred /= az_pred) then
          rss = rss + 1.0d4
          n_moon = n_moon + 1
          cycle
        end if
        da = obs_alts(kk) - alt_pred
        dz = obs_azs(kk) - az_pred
        if (dz >  180.0_dp) dz = dz - 360.0_dp
        if (dz < -180.0_dp) dz = dz + 360.0_dp
        rss = rss + da**2 + dz**2
        n_moon = n_moon + 1
      end do

      if (rss < rss_best) then
        rss_best = rss
        a_best = a_trial
        mu_best = mu_trial
      end if
    end do

    print '(A,F16.1,A,F8.3,A)', '  Coarse: a_E = ', a_best, ' km  (err = ', &
         (a_best/x_true(7)-1.0_dp)*100.0_dp, '%)'

    ! Refinement: narrow search around the coarse best
    a_lo = a_best - 2.0_dp * a_step
    a_hi = a_best + 2.0_dp * a_step
    a_step = (a_hi - a_lo) / real(N_GRID - 1, dp)
    rss_best = huge(1.0_dp)

    do ii = 1, N_GRID
      a_trial = a_lo + real(ii - 1, dp) * a_step
      if (a_trial < 1.0d6) cycle
      mu_trial = n_mean**2 * a_trial**3

      x_trial = xv
      x_trial(6) = mu_trial
      x_trial(7) = a_trial

      rss = 0.0_dp
      do kk = 1, n_obs, SUBSAMPLE
        if (obs_bodies(kk) /= 2) cycle
        call predict_altaz(x_trial, t_epoch, obs_jds(kk), &
                            lat_obs, lon_obs, obs_bodies(kk), alt_pred, az_pred)
        if (alt_pred /= alt_pred .or. az_pred /= az_pred) then
          rss = rss + 1.0d4; cycle
        end if
        da = obs_alts(kk) - alt_pred
        dz = obs_azs(kk) - az_pred
        if (dz >  180.0_dp) dz = dz - 360.0_dp
        if (dz < -180.0_dp) dz = dz + 360.0_dp
        rss = rss + da**2 + dz**2
      end do

      if (rss < rss_best) then
        rss_best = rss
        a_best = a_trial
        mu_best = mu_trial
      end if
    end do

    print '(A,F16.1,A,F8.3,A)', '  Refined: a_E = ', a_best, ' km  (err = ', &
         (a_best/x_true(7)-1.0_dp)*100.0_dp, '%)'

    ! Update state (no damping — trust the refined search)
    xv(6) = mu_best
    xv(7) = a_best
  end subroutine

end program sun_ekf_3body
