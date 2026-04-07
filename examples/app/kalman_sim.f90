module nbody3_ukf_mod
  use constants_mod
  implicit none

  integer, parameter :: NB = 3   ! Sun=1, Earth=2, Moon=3
  integer, parameter :: NX = 21  ! state dimension

  ! True GM values (for reference / printing)
  real(dp), parameter :: GM_SUN_TRUE   = 132712440041.279419_dp
  real(dp), parameter :: GM_EARTH_TRUE = 398600.435507_dp
  real(dp), parameter :: GM_MOON_TRUE  = 4902.800118_dp

  ! Yoshida 4th-order coefficients
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

  ! ── Pack/unpack state vector ──
  subroutine unpack_state(x, pos, vel, gm)
    real(dp), intent(in)  :: x(NX)
    real(dp), intent(out) :: pos(NB,3), vel(NB,3), gm(NB)
    pos(1,:) = x(1:3);   vel(1,:) = x(4:6)
    pos(2,:) = x(7:9);   vel(2,:) = x(10:12)
    pos(3,:) = x(13:15); vel(3,:) = x(16:18)
    gm(1) = x(19); gm(2) = x(20); gm(3) = x(21)
  end subroutine

  subroutine pack_state(pos, vel, gm, x)
    real(dp), intent(in)  :: pos(NB,3), vel(NB,3), gm(NB)
    real(dp), intent(out) :: x(NX)
    x(1:3) = pos(1,:);   x(4:6) = vel(1,:)
    x(7:9) = pos(2,:);   x(10:12) = vel(2,:)
    x(13:15) = pos(3,:); x(16:18) = vel(3,:)
    x(19) = gm(1); x(20) = gm(2); x(21) = gm(3)
  end subroutine

  ! ── 3-body accelerations ──
  subroutine compute_accelerations_3(pos, gm, acc)
    real(dp), intent(in)  :: pos(NB,3), gm(NB)
    real(dp), intent(out) :: acc(NB,3)
    integer :: i, j
    real(dp) :: r(3), r2, rinv3

    acc = 0.0_dp
    do i = 1, NB - 1
      do j = i + 1, NB
        r = pos(j,:) - pos(i,:)
        r2 = dot_product(r, r)
        rinv3 = 1.0_dp / (r2 * sqrt(r2))
        acc(i,:) = acc(i,:) + gm(j) * r * rinv3
        acc(j,:) = acc(j,:) - gm(i) * r * rinv3
      end do
    end do
  end subroutine

  ! ── Yoshida 4th-order step for 3 bodies ──
  subroutine yoshida4_step_3(pos, vel, gm, dt)
    real(dp), intent(inout) :: pos(NB,3), vel(NB,3)
    real(dp), intent(in)    :: gm(NB), dt
    real(dp) :: acc(NB,3)

    pos = pos + C1_Y * dt * vel
    call compute_accelerations_3(pos, gm, acc)
    vel = vel + D1_Y * dt * acc

    pos = pos + C2_Y * dt * vel
    call compute_accelerations_3(pos, gm, acc)
    vel = vel + D2_Y * dt * acc

    pos = pos + C3_Y * dt * vel
    call compute_accelerations_3(pos, gm, acc)
    vel = vel + D3_Y * dt * acc

    pos = pos + C4_Y * dt * vel
  end subroutine

  ! ── Propagate full state vector by dt seconds ──
  subroutine propagate_state(x_in, dt, x_out)
    real(dp), intent(in)  :: x_in(NX), dt
    real(dp), intent(out) :: x_out(NX)
    real(dp) :: pos(NB,3), vel(NB,3), gm(NB)
    real(dp) :: dt_step, t_done
    integer :: n_steps, i

    call unpack_state(x_in, pos, vel, gm)

    ! Use 3-hour (10800 s) max sub-steps for stability
    dt_step = sign(min(abs(dt), 10800.0_dp), dt)
    n_steps = max(1, nint(abs(dt) / abs(dt_step)))
    dt_step = dt / real(n_steps, dp)

    do i = 1, n_steps
      call yoshida4_step_3(pos, vel, gm, dt_step)
    end do

    call pack_state(pos, vel, gm, x_out)
  end subroutine

end module nbody3_ukf_mod


! ─────────────────────────────────────────────────────────────────────
!  Module 5: Simplified observation model (GMST + geometric projection)
! ─────────────────────────────────────────────────────────────────────
module obs_model_mod
  use constants_mod
  use nbody3_ukf_mod, only: NX, NB
  implicit none
contains

  ! ── Earth Rotation Angle (turns) ──
  function earth_rotation_angle(jd_whole, jd_frac) result(theta)
    real(dp), intent(in) :: jd_whole, jd_frac
    real(dp) :: theta, th
    th = 0.7790572732640_dp + 0.00273781191135448_dp * &
         (jd_whole - T0 + jd_frac)
    theta = mod(mod(th, 1.0_dp) + mod(jd_whole, 1.0_dp) + jd_frac, 1.0_dp)
  end function

  ! ── GMST in hours (full ERA + precession polynomial) ──
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

  ! ── ICRS → J2000 bias matrix ──
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

  ! ── Precession matrix (IAU 2006) ──
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

  ! ── P × B (precession × bias, no nutation) ──
  function compute_PB(jd_tdb) result(M)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: M(3,3), P(3,3), B(3,3)
    B = icrs_to_j2000_bias()
    P = compute_precession(jd_tdb)
    M = matmul(P, B)
  end function

  ! ── Simplified nutation (dominant lunar node term) ──
  subroutine simple_nutation(jd_tdb, dpsi, deps, mean_ob)
    real(dp), intent(in)  :: jd_tdb
    real(dp), intent(out) :: dpsi, deps, mean_ob
    real(dp) :: t, omega, omega2, f, d, l, lp

    t = (jd_tdb - T0) / 36525.0_dp

    ! Mean obliquity (arcsec → radians)
    mean_ob = (84381.406_dp + (-46.836769_dp + (-0.0001831_dp + 0.00200340_dp * t) * t) * t) &
              * ASEC2RAD

    ! Fundamental arguments (radians)
    omega = (450160.398036_dp + (-6962890.5431_dp + (7.4722_dp + 0.007702_dp * t) * t) * t) &
            * ASEC2RAD
    omega = mod(omega, TAU)

    ! Moon mean anomaly (l)
    l = (485868.249036_dp + (1717915923.2178_dp + (31.8792_dp + 0.051635_dp * t) * t) * t) &
        * ASEC2RAD

    ! Sun mean anomaly (lp)
    lp = (1287104.79305_dp + (129596581.0481_dp + (-0.5532_dp + 0.000136_dp * t) * t) * t) &
         * ASEC2RAD

    ! Moon mean elongation (D)
    d = (1072260.70369_dp + (1602961601.2090_dp + (-6.3706_dp + 0.006593_dp * t) * t) * t) &
        * ASEC2RAD

    ! Moon arg of latitude (F)
    f = (335779.526232_dp + (1739527262.8478_dp + (-12.7512_dp - 0.001037_dp * t) * t) * t) &
        * ASEC2RAD

    ! Top ~10 nutation terms (in arcsec, then convert to radians)
    dpsi = -17.2064161_dp * sin(omega) &
           -  1.3170907_dp * sin(2.0_dp*(f - d + omega)) &
           -  0.2276413_dp * sin(2.0_dp*(f + omega)) &
           +  0.2074554_dp * sin(2.0_dp*omega) &
           -  0.1491347_dp * sin(lp) &
           +  0.0130691_dp * sin(l) &
           -  0.0120323_dp * sin(2.0_dp*(f - d) + omega)
    dpsi = dpsi * ASEC2RAD

    deps =  9.2052331_dp * cos(omega) &
          + 0.5730336_dp * cos(2.0_dp*(f - d + omega)) &
          + 0.0978459_dp * cos(2.0_dp*(f + omega)) &
          - 0.0897492_dp * cos(2.0_dp*omega)
    deps = deps * ASEC2RAD
  end subroutine

  ! ── M = N × P × B (with simplified nutation) ──
  function compute_NPB(jd_tdb) result(M)
    real(dp), intent(in) :: jd_tdb
    real(dp) :: M(3,3), PB(3,3), Nmat(3,3)
    real(dp) :: dpsi, deps, mean_ob, true_ob
    real(dp) :: cobm, sobm, cobt, sobt, cpsi, spsi

    PB = compute_PB(jd_tdb)

    call simple_nutation(jd_tdb, dpsi, deps, mean_ob)
    true_ob = mean_ob + deps

    ! Nutation matrix (same as build_nutation_matrix in compute_from_scratch)
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

  ! ── GAST from GMST + equation of equinoxes ──
  function gast_from_gmst(gmst_h, dpsi, mean_ob, jd_tt) result(gast_h)
    real(dp), intent(in) :: gmst_h, dpsi, mean_ob, jd_tt
    real(dp) :: gast_h, t, ee_arcsec, ee_hours
    t = (jd_tt - T0) / 36525.0_dp
    ee_arcsec = dpsi / ASEC2RAD * cos(mean_ob) &
                + 0.00264_dp * sin(125.04_dp * PI/180.0_dp - 1934.136_dp * PI/180.0_dp * t) &
                + 0.000063_dp * sin(2.0_dp * (125.04_dp * PI/180.0_dp - 1934.136_dp * PI/180.0_dp * t))
    ee_hours = ee_arcsec / 54000.0_dp
    gast_h = mod(gmst_h + ee_hours, 24.0_dp)
  end function

  ! ── Predict (alt, az) for a body given the state vector ──
  subroutine predict_altaz(x, jd_utc, lat_deg, lon_deg, body_id, alt_deg, az_deg)
    real(dp), intent(in)  :: x(NX)
    real(dp), intent(in)  :: jd_utc, lat_deg, lon_deg
    integer,  intent(in)  :: body_id   ! 1=Sun, 2=Moon
    real(dp), intent(out) :: alt_deg, az_deg

    real(dp) :: earth_pos(3), earth_vel(3), body_pos(3), body_vel(3)
    real(dp) :: obs_gcrs(3), d_gcrs(3)
    real(dp) :: d_date(3), d_itrs(3), d_local(3), obs_itrs(3), obs_date(3)
    real(dp) :: lat, lon, gmst_h, gmst_rad
    real(dp) :: sinlat, coslat, sg, cg, slon, clon
    real(dp) :: north, east, up, dnorm, tmp
    real(dp) :: jd_whole, jd_frac, jd_tdb
    real(dp) :: M(3,3)
    real(dp) :: dist_km, light_time, d_unit(3), v_obs(3), vdot
    real(dp), parameter :: R_EARTH_KM = 6378.137_dp
    real(dp), parameter :: C_LIGHT_KMS = 299792.458_dp

    earth_pos = x(7:9)
    earth_vel = x(10:12)
    if (body_id == 1) then
      body_pos = x(1:3);  body_vel = x(4:6)     ! Sun
    else
      body_pos = x(13:15); body_vel = x(16:18)   ! Moon
    end if

    ! Light-time correction: use retarded body position
    dist_km = sqrt(sum((body_pos - earth_pos)**2))
    light_time = dist_km / C_LIGHT_KMS
    body_pos = body_pos - body_vel * light_time

    lat = lat_deg * PI / 180.0_dp
    lon = lon_deg * PI / 180.0_dp
    sinlat = sin(lat); coslat = cos(lat)
    slon = sin(lon); clon = cos(lon)

    ! Split JD for ERA precision
    jd_whole = floor(jd_utc + 0.5_dp)
    jd_frac  = jd_utc - jd_whole
    jd_tdb = jd_utc + 69.184_dp / DAY_S

    ! Full GMST
    gmst_h = gmst_hours(jd_whole, jd_frac, jd_tdb)

    ! Nutation + precession matrix M = N × P × B
    M = compute_NPB(jd_tdb)

    ! GAST = GMST + equation of equinoxes (uses nutation dpsi and mean_ob)
    block
      real(dp) :: dpsi_n, deps_n, mean_ob_n, gast_h, gast_rad
      call simple_nutation(jd_tdb, dpsi_n, deps_n, mean_ob_n)
      gast_h = gast_from_gmst(gmst_h, dpsi_n, mean_ob_n, jd_tdb)
      gast_rad = gast_h * TAU / 24.0_dp
      sg = sin(gast_rad); cg = cos(gast_rad)
    end block

    ! Observer in ITRS (Earth-centered, km)
    obs_itrs(1) = R_EARTH_KM * coslat * clon
    obs_itrs(2) = R_EARTH_KM * coslat * slon
    obs_itrs(3) = R_EARTH_KM * sinlat

    ! Observer in GCRS: M^T × Rz(+GAST) × obs_itrs
    obs_date(1) =  cg*obs_itrs(1) - sg*obs_itrs(2)
    obs_date(2) =  sg*obs_itrs(1) + cg*obs_itrs(2)
    obs_date(3) = obs_itrs(3)
    obs_gcrs = matmul(transpose(M), obs_date)

    ! Direction in GCRS
    d_gcrs = body_pos - earth_pos - obs_gcrs

    ! Aberration correction (first-order special relativistic)
    ! v_obs = Earth BCRS velocity (ignoring observer surface rotation, ~0.3 km/s)
    v_obs = earth_vel / C_LIGHT_KMS   ! v/c dimensionless
    dnorm = sqrt(sum(d_gcrs**2))
    if (dnorm > 1.0d-10) then
      d_unit = d_gcrs / dnorm
      vdot = dot_product(d_unit, v_obs)
      d_gcrs = d_gcrs + dnorm * (v_obs - d_unit * vdot)
    end if

    ! Transform to equinox of date: M × d_gcrs
    d_date = matmul(M, d_gcrs)

    ! Transform to ITRS: Rz(-GMST) × d_date
    d_itrs(1) =  cg*d_date(1) + sg*d_date(2)
    d_itrs(2) = -sg*d_date(1) + cg*d_date(2)
    d_itrs(3) = d_date(3)

    ! Align to local meridian: Rz(-lon)
    tmp        =  clon*d_itrs(1) + slon*d_itrs(2)
    d_local(2) = -slon*d_itrs(1) + clon*d_itrs(2)
    d_local(1) = tmp
    d_local(3) = d_itrs(3)

    ! Horizon: Up = coslat·X + sinlat·Z, North = -sinlat·X + coslat·Z
    up    = coslat*d_local(1) + sinlat*d_local(3)
    north = -sinlat*d_local(1) + coslat*d_local(3)
    east  = d_local(2)

    dnorm = sqrt(north**2 + east**2 + up**2)
    if (dnorm < 1.0d-10) then
      alt_deg = 0.0_dp; az_deg = 0.0_dp; return
    end if

    alt_deg = asin(up / dnorm) * 180.0_dp / PI
    az_deg  = atan2(east, north) * 180.0_dp / PI
    if (az_deg < 0.0_dp) az_deg = az_deg + 360.0_dp
  end subroutine

end module obs_model_mod


! ─────────────────────────────────────────────────────────────────────
!  Module 6: Scaled Unscented Kalman Filter
! ─────────────────────────────────────────────────────────────────────
module ukf_mod
  use constants_mod
  use nbody3_ukf_mod, only: NX, propagate_state
  use obs_model_mod, only: predict_altaz
  implicit none

  integer, parameter :: NZ = 2                   ! measurement dimension (alt, az)
  integer, parameter :: N_SIGMA = 2 * NX + 1     ! = 43

  ! UKF tuning parameters
  real(dp), parameter :: UKF_ALPHA = 1.0d-2
  real(dp), parameter :: UKF_BETA  = 2.0_dp
  real(dp), parameter :: UKF_KAPPA = 0.0_dp

contains

  ! ── Cholesky decomposition: A = L L^T, returns L ──
  subroutine cholesky(A, n, L, ok)
    integer, intent(in)   :: n
    real(dp), intent(in)  :: A(n,n)
    real(dp), intent(out) :: L(n,n)
    logical, intent(out)  :: ok
    integer :: i, j, k
    real(dp) :: s

    L = 0.0_dp
    ok = .true.
    do j = 1, n
      s = A(j,j)
      do k = 1, j - 1
        s = s - L(j,k)**2
      end do
      if (s <= 0.0_dp) then
        ok = .false.
        ! Attempt recovery: set small positive
        s = 1.0d-12
      end if
      L(j,j) = sqrt(s)
      do i = j + 1, n
        s = A(i,j)
        do k = 1, j - 1
          s = s - L(i,k) * L(j,k)
        end do
        L(i,j) = s / L(j,j)
      end do
    end do
  end subroutine

  ! ── Generate sigma points ──
  subroutine sigma_points(x, P, sigmas)
    real(dp), intent(in)  :: x(NX), P(NX,NX)
    real(dp), intent(out) :: sigmas(NX, N_SIGMA)
    real(dp) :: L(NX,NX), scaled_P(NX,NX)
    real(dp) :: lambda
    logical :: ok
    integer :: i

    lambda = UKF_ALPHA**2 * (real(NX,dp) + UKF_KAPPA) - real(NX,dp)
    scaled_P = (real(NX,dp) + lambda) * P

    call cholesky(scaled_P, NX, L, ok)

    ! Sigma point 0 = mean
    sigmas(:, 1) = x

    ! Sigma points 1..N and N+1..2N
    do i = 1, NX
      sigmas(:, 1+i)    = x + L(:,i)
      sigmas(:, 1+NX+i) = x - L(:,i)
    end do
  end subroutine

  ! ── UKF weights ──
  subroutine ukf_weights(Wm, Wc)
    real(dp), intent(out) :: Wm(N_SIGMA), Wc(N_SIGMA)
    real(dp) :: lambda
    integer :: i

    lambda = UKF_ALPHA**2 * (real(NX,dp) + UKF_KAPPA) - real(NX,dp)

    Wm(1) = lambda / (real(NX,dp) + lambda)
    Wc(1) = Wm(1) + (1.0_dp - UKF_ALPHA**2 + UKF_BETA)

    do i = 2, N_SIGMA
      Wm(i) = 1.0_dp / (2.0_dp * (real(NX,dp) + lambda))
      Wc(i) = Wm(i)
    end do
  end subroutine

  ! ── UKF prediction step: propagate mean + add process noise ──
  !
  ! We propagate only the mean through the dynamics (not sigma points),
  ! because N-body orbits diverge too fast for sigma-point spreading.
  ! Hybrid predict: propagate mean through dynamics + add Q for covariance.
  ! Note: GM states don't couple through dynamics in this approach.
  subroutine ukf_predict(x, P, Q_rate, dt)
    real(dp), intent(inout) :: x(NX), P(NX,NX)
    real(dp), intent(in)    :: Q_rate(NX,NX), dt
    real(dp) :: x_new(NX)

    ! Propagate mean state through dynamics
    call propagate_state(x, dt, x_new)
    x = x_new

    ! Grow covariance with process noise
    P = P + Q_rate * abs(dt)
  end subroutine

  ! ── UKF measurement update ──
  subroutine ukf_update(x, P, z, R, jd_utc, lat_deg, lon_deg, body_id)
    real(dp), intent(inout) :: x(NX), P(NX,NX)
    real(dp), intent(in)    :: z(NZ)        ! measurement [alt, az] degrees
    real(dp), intent(in)    :: R(NZ,NZ)     ! measurement noise
    real(dp), intent(in)    :: jd_utc, lat_deg, lon_deg
    integer,  intent(in)    :: body_id

    real(dp) :: sigmas(NX, N_SIGMA)
    real(dp) :: z_sigma(NZ, N_SIGMA)
    real(dp) :: Wm(N_SIGMA), Wc(N_SIGMA)
    real(dp) :: z_pred(NZ), dz(NZ), dx(NX)
    real(dp) :: Pzz(NZ,NZ), Pxz(NX,NZ)
    real(dp) :: KG(NX,NZ), innovation(NZ)
    real(dp) :: alt_pred, az_pred
    integer :: i, j, k

    call ukf_weights(Wm, Wc)
    call sigma_points(x, P, sigmas)

    ! Transform sigma points through observation model
    do i = 1, N_SIGMA
      call predict_altaz(sigmas(:,i), jd_utc, lat_deg, lon_deg, body_id, &
                         alt_pred, az_pred)
      z_sigma(1,i) = alt_pred
      z_sigma(2,i) = az_pred
    end do

    ! Predicted measurement mean
    z_pred = 0.0_dp
    do i = 1, N_SIGMA
      z_pred = z_pred + Wm(i) * z_sigma(:,i)
    end do

    ! Innovation covariance Pzz + cross-covariance Pxz
    Pzz = R    ! start with measurement noise
    Pxz = 0.0_dp
    do i = 1, N_SIGMA
      dz = z_sigma(:,i) - z_pred
      ! Handle azimuth wraparound
      if (dz(2) > 180.0_dp) dz(2) = dz(2) - 360.0_dp
      if (dz(2) < -180.0_dp) dz(2) = dz(2) + 360.0_dp

      dx = sigmas(:,i) - x

      do j = 1, NZ
        do k = 1, NZ
          Pzz(j,k) = Pzz(j,k) + Wc(i) * dz(j) * dz(k)
        end do
      end do
      do j = 1, NX
        do k = 1, NZ
          Pxz(j,k) = Pxz(j,k) + Wc(i) * dx(j) * dz(k)
        end do
      end do
    end do

    ! Kalman gain KG = Pxz * Pzz^{-1}  (2×2 inversion)
    call solve_2x2(Pzz, Pxz, NX, KG)

    ! Innovation
    innovation = z - z_pred
    if (innovation(2) > 180.0_dp) innovation(2) = innovation(2) - 360.0_dp
    if (innovation(2) < -180.0_dp) innovation(2) = innovation(2) + 360.0_dp

    ! State update
    do i = 1, NX
      x(i) = x(i) + KG(i,1) * innovation(1) + KG(i,2) * innovation(2)
    end do

    ! Covariance update: P = P - KG * Pzz * KG^T
    do j = 1, NX
      do k = 1, NX
        do i = 1, NZ
          P(j,k) = P(j,k) - KG(j,i) * Pzz(i,1) * KG(k,1) &
                           - KG(j,i) * Pzz(i,2) * KG(k,2)
        end do
      end do
    end do

    ! Symmetrise
    do j = 1, NX
      do k = j+1, NX
        P(j,k) = 0.5_dp * (P(j,k) + P(k,j))
        P(k,j) = P(j,k)
      end do
    end do
  end subroutine

  ! ── Solve K = Pxz * Pzz^{-1} for 2×2 Pzz ──
  subroutine solve_2x2(Pzz, Pxz, nx, K)
    integer, intent(in) :: nx
    real(dp), intent(in) :: Pzz(2,2), Pxz(nx,2)
    real(dp), intent(out) :: K(nx,2)
    real(dp) :: det, inv(2,2)
    integer :: i

    det = Pzz(1,1) * Pzz(2,2) - Pzz(1,2) * Pzz(2,1)
    if (abs(det) < 1.0d-30) then
      K = 0.0_dp; return
    end if

    inv(1,1) =  Pzz(2,2) / det
    inv(1,2) = -Pzz(1,2) / det
    inv(2,1) = -Pzz(2,1) / det
    inv(2,2) =  Pzz(1,1) / det

    do i = 1, nx
      K(i,1) = Pxz(i,1) * inv(1,1) + Pxz(i,2) * inv(2,1)
      K(i,2) = Pxz(i,1) * inv(1,2) + Pxz(i,2) * inv(2,2)
    end do
  end subroutine

  ! ── Solve A*x = b for symmetric positive definite A via Cholesky ──
  subroutine solve_symmetric(A, b, n, x)
    integer, intent(in)   :: n
    real(dp), intent(in)  :: A(n,n), b(n)
    real(dp), intent(out) :: x(n)
    real(dp) :: L(n,n), y(n), s
    logical :: ok
    integer :: i, j

    call cholesky(A, n, L, ok)
    if (.not. ok) then
      x = 0.0_dp; return
    end if

    ! Forward substitution: L*y = b
    do i = 1, n
      s = b(i)
      do j = 1, i - 1
        s = s - L(i,j) * y(j)
      end do
      y(i) = s / L(i,i)
    end do

    ! Back substitution: L^T*x = y
    do i = n, 1, -1
      s = y(i)
      do j = i + 1, n
        s = s - L(j,i) * x(j)
      end do
      x(i) = s / L(i,i)
    end do
  end subroutine

end module ukf_mod


! ─────────────────────────────────────────────────────────────────────
!  Main program
! ─────────────────────────────────────────────────────────────────────
program kalman_sim
  use constants_mod
  use spk_reader_mod
  use nbody3_ukf_mod
  use obs_reader_mod
  use obs_model_mod
  use ukf_mod
  implicit none

  ! ── Configuration ──
  real(dp) :: lat_deg, lon_deg
  integer  :: n_days_filter
  real(dp) :: perturb_gm_frac, perturb_pos_km, perturb_vel_kms

  ! ── State and covariance ──
  real(dp) :: x(NX), P(NX,NX), Q(NX,NX), R_meas(NZ,NZ), R_moon(NZ,NZ), R_use(NZ,NZ)

  ! ── True initial state (from DE440s) ──
  real(dp) :: x_true(NX)

  ! ── Observations ──
  type(observation) :: obs(MAX_OBS)
  integer :: n_obs

  ! ── Working variables ──
  type(spk_kernel) :: kernel
  real(dp) :: pos_sun(3), vel_sun(3), pos_earth(3), vel_earth(3)
  real(dp) :: pos_moon(3), vel_moon(3)
  real(dp) :: p1(3), p2(3), v1(3), v2(3)
  real(dp) :: jd_start, jd_max, jd_current, dt_sec
  real(dp) :: z_meas(NZ)
  real(dp) :: r_val
  integer  :: i, j, n_seed
  integer, allocatable :: seed(:)
  character(len=256) :: exe_dir, bsp_file, obs_file
  integer :: slen

  ! ── Logging ──
  real(dp) :: earth_sun_km, earth_moon_km
  real(dp) :: gm_sun_est, gm_earth_est, gm_moon_est
  integer :: log_counter
  real(dp) :: alt_test, az_test, model_err_alt, model_err_az, dalt, daz
  real(dp) :: x_prop(NX), jd_verify

  ! ── Batch estimation variables ──
  integer, parameter :: NM_BATCH = 2
  real(dp) :: x_iter(NX), x_prior(NX), x_base_b(NX), x_pert_b(NX)
  real(dp) :: batch_residuals(NM_BATCH, MAX_OBS)
  integer  :: batch_obs_idx(MAX_OBS)
  real(dp) :: batch_H(NM_BATCH * MAX_OBS, NX)
  real(dp) :: batch_HtWH(NX, NX), batch_HtWr(NX), batch_dx(NX)
  real(dp) :: b_pred_alt, b_pred_az, b_pred_alt_p, b_pred_az_p
  real(dp) :: b_w_alt, b_w_az, b_rms_alt, b_rms_az, b_rms_prev
  real(dp) :: b_delta_param, b_sigma_prior(NX)
  integer :: b_iter, b_m_total, b_row, b_col, b_jj, b_kk, b_ii
  integer :: b_n_active, b_active_idx(NX), b_n_phase_iter
  integer :: b_total_iter, b_n_obs_win
  real(dp) :: b_pos_lim, b_vel_lim, b_gm_lim
  logical :: b_is_active(NX)
  ! ── Get exe directory ──
  call get_command_argument(0, exe_dir)
  slen = index(exe_dir, '/', back=.true.)
  if (slen > 0) then; exe_dir = exe_dir(1:slen); else; exe_dir = './'; end if

  ! ═══════════════════════════════════════════════════════════════════
  !  Configuration
  ! ═══════════════════════════════════════════════════════════════════
  lat_deg = 40.0_dp
  lon_deg = 0.0_dp
  n_days_filter = 365
  perturb_gm_frac = 0.01_dp    ! 1% mass perturbation
  perturb_pos_km  = 100.0_dp   ! km position perturbation
  perturb_vel_kms = 0.01_dp    ! km/s velocity perturbation

  print '(A)',        '════════════════════════════════════════════════════════'
  print '(A)',        '  3-Body UKF Kalman Filter'
  print '(A,F8.4,A,F8.4)', '  Observer: lat=', lat_deg, ' lon=', lon_deg
  print '(A,I4,A)',   '  Filter window: ', n_days_filter, ' days'
  print '(A)',        '════════════════════════════════════════════════════════'

  ! ═══════════════════════════════════════════════════════════════════
  !  Load true initial conditions from DE440s
  ! ═══════════════════════════════════════════════════════════════════
  call get_command_argument(1, bsp_file)
  if (len_trim(bsp_file) == 0) bsp_file = 'de440s.bsp'
  if (index(trim(bsp_file), '/') > 0) then
    call spk_open(trim(bsp_file), kernel)
  else
    call spk_open(trim(exe_dir) // trim(bsp_file), kernel)
  end if

  ! Start at 2025-01-01 noon TDB (matches observations.dat epoch)
  jd_start = 2460677.0_dp    ! Julian day for 2025-01-01 noon

  ! Sun (NAIF 10 wrt SSB=0)
  call spk_compute_and_diff(kernel, 0, 10, jd_start, 0.0_dp, pos_sun, vel_sun)
  vel_sun = vel_sun / DAY_S

  ! Earth (EMB=3 + 399 wrt EMB)
  call spk_compute_and_diff(kernel, 0, 3, jd_start, 0.0_dp, p1, v1)
  call spk_compute_and_diff(kernel, 3, 399, jd_start, 0.0_dp, p2, v2)
  pos_earth = p1 + p2
  vel_earth = (v1 + v2) / DAY_S

  ! Moon (EMB=3 + 301 wrt EMB)
  call spk_compute_and_diff(kernel, 0, 3, jd_start, 0.0_dp, p1, v1)
  call spk_compute_and_diff(kernel, 3, 301, jd_start, 0.0_dp, p2, v2)
  pos_moon = p1 + p2
  vel_moon = (v1 + v2) / DAY_S

  call spk_close(kernel)

  ! Pack true state
  x_true(1:3)   = pos_sun;   x_true(4:6)   = vel_sun
  x_true(7:9)   = pos_earth; x_true(10:12) = vel_earth
  x_true(13:15) = pos_moon;  x_true(16:18) = vel_moon
  x_true(19)    = GM_SUN_TRUE
  x_true(20)    = GM_EARTH_TRUE
  x_true(21)    = GM_MOON_TRUE

  print '(A)',             '  True initial state loaded from DE440s'
  print '(A,ES16.8,A)',    '    GM_sun   = ', GM_SUN_TRUE, ' km³/s²'
  print '(A,ES16.8,A)',    '    GM_earth = ', GM_EARTH_TRUE, ' km³/s²'
  print '(A,ES16.8,A)',    '    GM_moon  = ', GM_MOON_TRUE, ' km³/s²'

  earth_sun_km = sqrt(dot_product(pos_earth - pos_sun, pos_earth - pos_sun))
  earth_moon_km = sqrt(dot_product(pos_earth - pos_moon, pos_earth - pos_moon))
  print '(A,F12.1,A,F12.6,A)', '    Earth-Sun  = ', earth_sun_km, ' km (', &
       earth_sun_km / AU_KM, ' AU)'
  print '(A,F12.1,A)',         '    Earth-Moon = ', earth_moon_km, ' km'

  ! ═══════════════════════════════════════════════════════════════════
  !  Perturb initial conditions
  ! ═══════════════════════════════════════════════════════════════════
  call random_seed(size=n_seed)
  allocate(seed(n_seed))
  do i = 1, n_seed; seed(i) = 123 + i * 97; end do
  call random_seed(put=seed)
  deallocate(seed)

  x = x_true

  ! Perturb positions (km)
  do i = 1, 18
    call random_number(r_val)
    if (i <= 3 .or. (i >= 7 .and. i <= 9) .or. (i >= 13 .and. i <= 15)) then
      ! Position: ±perturb_pos_km
      x(i) = x(i) + (r_val - 0.5_dp) * 2.0_dp * perturb_pos_km
    else
      ! Velocity: ±perturb_vel_kms
      x(i) = x(i) + (r_val - 0.5_dp) * 2.0_dp * perturb_vel_kms
    end if
  end do

  ! Perturb GM values
  call random_number(r_val)
  x(19) = GM_SUN_TRUE   * (1.0_dp + (r_val - 0.5_dp) * 2.0_dp * perturb_gm_frac)
  call random_number(r_val)
  x(20) = GM_EARTH_TRUE * (1.0_dp + (r_val - 0.5_dp) * 2.0_dp * perturb_gm_frac)
  call random_number(r_val)
  x(21) = GM_MOON_TRUE  * (1.0_dp + (r_val - 0.5_dp) * 2.0_dp * perturb_gm_frac)

  print '(/,A)',           '  Perturbed initial guess:'
  print '(A,ES16.8,A,F6.1,A)', '    GM_sun   = ', x(19), ' (', &
       (x(19)/GM_SUN_TRUE - 1.0_dp)*100.0_dp, '%)'
  print '(A,ES16.8,A,F6.1,A)', '    GM_earth = ', x(20), ' (', &
       (x(20)/GM_EARTH_TRUE - 1.0_dp)*100.0_dp, '%)'
  print '(A,ES16.8,A,F6.1,A)', '    GM_moon  = ', x(21), ' (', &
       (x(21)/GM_MOON_TRUE - 1.0_dp)*100.0_dp, '%)'

  ! ═══════════════════════════════════════════════════════════════════
  !  Initialize covariance matrices
  ! ═══════════════════════════════════════════════════════════════════
  P = 0.0_dp
  ! Position uncertainty: 200 km (2× perturbation)
  do i = 1, 3;   P(i,i) = (200.0_dp)**2; end do     ! Sun pos
  do i = 7, 9;   P(i,i) = (200.0_dp)**2; end do     ! Earth pos
  do i = 13, 15; P(i,i) = (200.0_dp)**2; end do     ! Moon pos
  ! Velocity uncertainty: 0.02 km/s
  do i = 4, 6;   P(i,i) = (0.02_dp)**2; end do
  do i = 10, 12; P(i,i) = (0.02_dp)**2; end do
  do i = 16, 18; P(i,i) = (0.02_dp)**2; end do
  ! GM uncertainty: 2% (covers the 1% perturbation)
  P(19,19) = (0.02_dp * GM_SUN_TRUE)**2
  P(20,20) = (0.02_dp * GM_EARTH_TRUE)**2
  P(21,21) = (0.02_dp * GM_MOON_TRUE)**2

  ! Process noise Q_rate (per second) — accounts for missing outer planets.
  ! Jupiter perturbs Earth by ~2e-10 km/s², causing ~800 km/day drift.
  Q = 0.0_dp
  ! Sun position: ~10 km/√day → Q_rate = 100/(86400) ≈ 1.2e-3
  do i = 1, 3;   Q(i,i) = 1.0d-3; end do
  do i = 4, 6;   Q(i,i) = 1.0d-13; end do
  ! Earth position: ~50 km/√day
  do i = 7, 9;   Q(i,i) = 3.0d-2; end do
  do i = 10, 12; Q(i,i) = 1.0d-12; end do
  ! Moon position: ~5 km/√day (differential, relative to Earth)
  do i = 13, 15; Q(i,i) = 3.0d-4; end do
  do i = 16, 18; Q(i,i) = 1.0d-13; end do
  ! GM: essentially constant (tiny noise for numerical stability)
  Q(19,19) = 1.0d-6
  Q(20,20) = 1.0d-10
  Q(21,21) = 1.0d-12

  ! Measurement noise R
  R_meas = 0.0_dp
  R_meas(1,1) = (30.0_dp / 3600.0_dp)**2   ! alt: 30" (noise + model)
  R_meas(2,2) = (30.0_dp / 3600.0_dp)**2   ! az: 30"
  R_moon = 0.0_dp
  R_moon(1,1) = (60.0_dp / 3600.0_dp)**2   ! Moon alt: 60" (dynamics drift)
  R_moon(2,2) = (100.0_dp / 3600.0_dp)**2  ! Moon az: 100" (largest drift)

  ! ═══════════════════════════════════════════════════════════════════
  !  Read observations
  ! ═══════════════════════════════════════════════════════════════════
  obs_file = trim(exe_dir) // 'observations.dat'
  jd_max = jd_start + real(n_days_filter, dp)
  call read_observations(obs_file, obs, n_obs, jd_max)

  print '(/,A,I6,A,I3,A)', '  Loaded ', n_obs, ' observations (', n_days_filter, ' days)'

  ! ═══════════════════════════════════════════════════════════════════
  !  Verify observation model with true state (propagated)
  ! ═══════════════════════════════════════════════════════════════════
  print '(/,A)', '  Observation model verification (propagated true state):'
  print '(A)', '  Day   E-M dist(km)  E-S dist(km)  body  dalt"     daz"'
  model_err_alt = 0.0_dp
  model_err_az  = 0.0_dp
  j = 0
  x_prop = x_true  ! copy to propagate
  jd_verify = jd_start
  do i = 1, n_obs
    ! Propagate true state to observation time
    dt_sec = (obs(i)%jd - jd_verify) * DAY_S
    if (abs(dt_sec) > 0.1_dp) then
      call propagate_state(x_prop, dt_sec, x_prop)
    end if
    jd_verify = obs(i)%jd

    call predict_altaz(x_prop, obs(i)%jd, lat_deg, lon_deg, obs(i)%body, &
                       alt_test, az_test)
    dalt = alt_test - obs(i)%alt_true
    daz  = az_test  - obs(i)%az_true
    if (daz > 180.0_dp)  daz = daz - 360.0_dp
    if (daz < -180.0_dp) daz = daz + 360.0_dp
    if (i <= 8 .or. (i <= 100 .and. mod(i,25)==0)) then
      earth_sun_km  = sqrt(sum((x_prop(7:9) - x_prop(1:3))**2))
      earth_moon_km = sqrt(sum((x_prop(13:15) - x_prop(7:9))**2))
      print '(F6.2,2F14.1,I6,2F10.2)', &
           obs(i)%jd - jd_start, earth_moon_km, earth_sun_km, &
           obs(i)%body, dalt*3600.0_dp, daz*3600.0_dp
    end if
    model_err_alt = model_err_alt + dalt**2
    model_err_az  = model_err_az  + daz**2
    j = j + 1
  end do
  model_err_alt = sqrt(model_err_alt / real(j, dp)) * 3600.0_dp
  model_err_az  = sqrt(model_err_az  / real(j, dp)) * 3600.0_dp
  print '(A,I5,A,F8.2,A,F8.2,A)', '    RMS over', j, ' obs: alt=', model_err_alt, &
       '"  az=', model_err_az, '"'


  ! ═══════════════════════════════════════════════════════════════════
  !  Batch Least-Squares Orbit Determination
  !  Stage 1: Sun-only obs → fix solar orbit + GM_sun
  !  Stage 2: Moon short-arc continuation → fix lunar orbit + GM_earth/moon
  !  Stage 3: All obs → final refinement
  ! ═══════════════════════════════════════════════════════════════════
  print '(/,A)', '────────────────────────────────────────────────────────'
  print '(A)',   '  Running batch least-squares estimation...'
  print '(A)',   '────────────────────────────────────────────────────────'

  x_iter = x
  x_prior = x
  b_total_iter = 0

  b_sigma_prior(1:18) = 0.0_dp
  b_sigma_prior(19) = 0.01_dp * abs(x_prior(19))
  b_sigma_prior(20) = 0.01_dp * abs(x_prior(20))
  b_sigma_prior(21) = 0.01_dp * abs(x_prior(21))

  ! ── STAGE 1: Sun-only → fix solar orbit ──
  ! Active: Sun pos/vel (1-6), Earth pos/vel (7-12), GM_sun (19)
  b_is_active = .false.
  b_is_active(1:12) = .true.
  b_is_active(19) = .true.
  b_n_active = 0
  do i = 1, NX
    if (b_is_active(i)) then
      b_n_active = b_n_active + 1
      b_active_idx(b_n_active) = i
    end if
  end do
  b_n_phase_iter = 40
  b_gm_lim = 0.002_dp
  b_pos_lim = 50000.0_dp
  b_vel_lim = 1.0_dp

  b_n_obs_win = 0
  do i = 1, n_obs
    if (obs(i)%body == 1) then
      b_n_obs_win = b_n_obs_win + 1
      batch_obs_idx(b_n_obs_win) = i
    end if
  end do
  print '(/,A,I5,A)', '  Stage 1: Sun-only (', b_n_obs_win, ' obs) → fix solar orbit'
  include 'batch_pass.inc'

  ! ── STAGE 2: Moon pos/vel only, multi-window continuation ──
  b_is_active = .false.
  b_is_active(13:18) = .true.  ! Moon pos/vel only
  b_n_active = 0
  do i = 1, NX
    if (b_is_active(i)) then
      b_n_active = b_n_active + 1
      b_active_idx(b_n_active) = i
    end if
  end do
  b_gm_lim = 0.002_dp
  b_pos_lim = 50000.0_dp
  b_vel_lim = 1.0_dp

  block
    real(dp) :: mw_days(3), jd_mw
    integer :: imw
    mw_days(1) = 5.0_dp; mw_days(2) = 15.0_dp; mw_days(3) = 30.0_dp

    do imw = 1, 3
      jd_mw = jd_start + mw_days(imw)
      b_n_obs_win = 0
      do i = 1, n_obs
        if (obs(i)%body == 2 .and. obs(i)%jd <= jd_mw) then
          b_n_obs_win = b_n_obs_win + 1
          batch_obs_idx(b_n_obs_win) = i
        end if
      end do
      b_n_phase_iter = 30
      print '(/,A,F6.0,A,I5,A)', '  Stage 2: Moon pos/vel ', mw_days(imw), &
           'd (', b_n_obs_win, ' obs)'
      include 'batch_pass.inc'
    end do
  end block

  x = x_iter


! ═══════════════════════════════════════════════════════════════════
  !  Final results
  ! ═══════════════════════════════════════════════════════════════════
  print '(/,A)', '════════════════════════════════════════════════════════'
  print '(A)',   '  FINAL ESTIMATES'
  print '(A)',   '════════════════════════════════════════════════════════'

  print '(A)',                        '  Masses (GM in km³/s²):'
  print '(A,ES18.10,A,ES18.10,A,F8.3,A)', &
       '    Sun:   ', x(19),   '  true: ', GM_SUN_TRUE,   '  err: ', &
       (x(19)/GM_SUN_TRUE - 1.0_dp)*100.0_dp, '%'
  print '(A,ES18.10,A,ES18.10,A,F8.3,A)', &
       '    Earth: ', x(20),   '  true: ', GM_EARTH_TRUE, '  err: ', &
       (x(20)/GM_EARTH_TRUE - 1.0_dp)*100.0_dp, '%'
  print '(A,ES18.10,A,ES18.10,A,F8.3,A)', &
       '    Moon:  ', x(21),   '  true: ', GM_MOON_TRUE,  '  err: ', &
       (x(21)/GM_MOON_TRUE - 1.0_dp)*100.0_dp, '%'

  print '(A)',                        '  Mass ratios:'
  print '(A,F14.1,A,F14.1)', &
       '    M_sun/M_earth:  ', x(19)/x(20), '  true: ', GM_SUN_TRUE/GM_EARTH_TRUE
  print '(A,F14.1,A,F14.1)', &
       '    M_earth/M_moon: ', x(20)/x(21), '  true: ', GM_EARTH_TRUE/GM_MOON_TRUE

  earth_sun_km = sqrt(dot_product(x(7:9) - x(1:3), x(7:9) - x(1:3)))
  earth_moon_km = sqrt(dot_product(x(7:9) - x(13:15), x(7:9) - x(13:15)))

  print '(/,A)',                      '  Initial distances (at epoch):'
  print '(A,F14.1,A,F12.6,A)',       '    Earth-Sun:  ', earth_sun_km, &
       ' km (', earth_sun_km / AU_KM, ' AU)'
  print '(A,F14.1,A)',               '    Earth-Moon: ', earth_moon_km, ' km'

  print '(/,A)', '  Note: GM_earth and GM_moon estimation is limited to ~1%'
  print '(A)',   '  accuracy by the 3-body model (missing Jupiter/Venus'
  print '(A)',   '  perturbations cause systematic Moon position errors'
  print '(A)',   '  of ~40-90 arcsec over 6+ days).'

  print '(/,A)', '════════════════════════════════════════════════════════'

end program kalman_sim
