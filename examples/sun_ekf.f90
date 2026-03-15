! ═══════════════════════════════════════════════════════════════════════
! sun_ekf.f90 — 2-body Sun-Earth EKF with N-body propagation
! ═══════════════════════════════════════════════════════════════════════
!
! State vector: x = [e, i, Ω, ω, M₀, μ, a]  (7 elements, all constant)
!   e    = eccentricity
!   i    = inclination (≈ obliquity ≈ 23.44° in equatorial J2000)
!   Ω    = longitude of ascending node
!   ω    = argument of periapsis
!   M₀   = mean anomaly at reference epoch
!   μ    = GM_sun + GM_earth (total gravitational parameter)
!   a    = semi-major axis (km)
!
! Fixed parameters:
!   R_E  = Earth radius (6378.137 km)
!
! Dynamics: F = I (elements are constants of motion in 2-body)
! Propagation: Yoshida 4th-order N-body integrator (N=2)
!   Elements → Cartesian at epoch → propagate → alt/az
! Observations: Sun alt/az from observations.dat
! Filter: EKF with finite-difference Jacobian
!
! ═══════════════════════════════════════════════════════════════════════

! ─── Module 1: Constants ─────────────────────────────────────────────
module constants_mod
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), parameter :: PI  = 3.14159265358979323846_dp
  real(dp), parameter :: TAU = 2.0_dp * PI
  real(dp), parameter :: DEG2RAD  = PI / 180.0_dp
  real(dp), parameter :: RAD2DEG  = 180.0_dp / PI
  real(dp), parameter :: ASEC2RAD = PI / (180.0_dp * 3600.0_dp)
  real(dp), parameter :: T0    = 2451545.0_dp      ! J2000.0
  real(dp), parameter :: DAY_S = 86400.0_dp
  real(dp), parameter :: AU_KM = 149597870.7_dp
  real(dp), parameter :: C_LIGHT_KMS = 299792.458_dp
  real(dp), parameter :: R_EARTH_KM  = 6378.137_dp
end module constants_mod

! ─── Module 2: SPK reader (DE440s) ──────────────────────────────────
module spk_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_SEGMENTS = 64

  type :: spk_segment
    real(dp) :: start_second, end_second
    integer  :: target, center, frame, data_type
    integer  :: start_i, end_i
    real(dp) :: start_jd, end_jd
    logical  :: loaded = .false.
    real(dp) :: init_epoch, intlen
    integer  :: n_intervals, coefficient_count, component_count
    real(dp), allocatable :: coeffs(:,:,:)
  end type

  type :: spk_kernel
    integer :: unit_num = -1
    integer :: n_segments = 0
    type(spk_segment) :: segments(MAX_SEGMENTS)
  end type

contains

  subroutine spk_open(filename, kernel)
    character(len=*), intent(in) :: filename
    type(spk_kernel), intent(out) :: kernel
    integer :: u, fward, nd, ni
    character(8) :: locidw
    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='old', action='read')
    kernel%unit_num = u
    read(u) locidw
    read(u) nd
    read(u) ni
    read(u, pos=77) fward
    call parse_summaries(u, fward, nd, ni, kernel)
  end subroutine

  subroutine parse_summaries(u, fward, nd, ni, kernel)
    integer, intent(in) :: u, fward, nd, ni
    type(spk_kernel), intent(inout) :: kernel
    integer :: record_number, n_summaries, i, seg_idx
    real(dp) :: next_rec, prev_rec, nsumm_d
    real(dp) :: start_sec, end_sec
    integer :: tgt, ctr, frm, dtype, si, ei
    integer :: base_pos, summary_size, step, ctrl_size

    summary_size = nd * 8 + ni * 4
    step = summary_size
    if (mod(step, 8) /= 0) step = step + (8 - mod(step, 8))
    ctrl_size = 24
    record_number = fward
    seg_idx = 0

    do while (record_number /= 0)
      base_pos = (record_number - 1) * 1024 + 1
      read(u, pos=base_pos) next_rec, prev_rec, nsumm_d
      n_summaries = int(nsumm_d)
      do i = 0, n_summaries - 1
        seg_idx = seg_idx + 1
        if (seg_idx > MAX_SEGMENTS) stop 'Too many segments'
        read(u, pos=base_pos + ctrl_size + i * step) start_sec, end_sec, &
             tgt, ctr, frm, dtype, si, ei
        kernel%segments(seg_idx)%start_second = start_sec
        kernel%segments(seg_idx)%end_second   = end_sec
        kernel%segments(seg_idx)%target       = tgt
        kernel%segments(seg_idx)%center       = ctr
        kernel%segments(seg_idx)%frame        = frm
        kernel%segments(seg_idx)%data_type    = dtype
        kernel%segments(seg_idx)%start_i      = si
        kernel%segments(seg_idx)%end_i        = ei
        kernel%segments(seg_idx)%start_jd     = T0 + start_sec / DAY_S
        kernel%segments(seg_idx)%end_jd       = T0 + end_sec / DAY_S
      end do
      record_number = int(next_rec)
    end do
    kernel%n_segments = seg_idx
  end subroutine

  function find_segment(kernel, center, target) result(idx)
    type(spk_kernel), intent(in) :: kernel
    integer, intent(in) :: center, target
    integer :: idx, i
    idx = -1
    do i = kernel%n_segments, 1, -1
      if (kernel%segments(i)%center == center .and. &
          kernel%segments(i)%target == target) then
        idx = i; return
      end if
    end do
  end function

  subroutine load_segment_data(kernel, idx)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: idx
    type(spk_segment) :: seg
    real(dp) :: meta(4)
    integer :: rsize_i, n_i, coeff_count, comp_count
    integer :: u, pos, total_words
    real(dp), allocatable :: raw(:)
    integer :: rec, c, k, kk

    seg = kernel%segments(idx)
    u = kernel%unit_num
    if (seg%loaded) return

    pos = (seg%end_i - 4) * 8 + 1
    read(u, pos=pos) meta
    kernel%segments(idx)%init_epoch = meta(1)
    kernel%segments(idx)%intlen     = meta(2)
    rsize_i = int(meta(3))
    n_i     = int(meta(4))
    kernel%segments(idx)%n_intervals = n_i

    if (seg%data_type == 2) then; comp_count = 3; else; comp_count = 6; end if
    coeff_count = (rsize_i - 2) / comp_count
    kernel%segments(idx)%coefficient_count = coeff_count
    kernel%segments(idx)%component_count   = comp_count

    total_words = rsize_i * n_i
    allocate(raw(total_words))
    pos = (seg%start_i - 1) * 8 + 1
    read(u, pos=pos) raw

    allocate(kernel%segments(idx)%coeffs(coeff_count, comp_count, n_i))
    do rec = 1, n_i
      do c = 1, comp_count
        do k = 1, coeff_count
          kk = coeff_count - k + 1
          kernel%segments(idx)%coeffs(kk, c, rec) = &
              raw((rec-1)*rsize_i + 2 + (c-1)*coeff_count + k)
        end do
      end do
    end do
    deallocate(raw)
    kernel%segments(idx)%loaded = .true.
  end subroutine

  subroutine spk_compute_and_diff(kernel, center, target, tdb_whole, tdb_frac, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: center, target
    real(dp), intent(in) :: tdb_whole, tdb_frac
    real(dp), intent(out) :: pos(3), vel(3)
    integer :: idx, n_int, cc
    real(dp) :: init_e, intlen
    real(dp) :: index1, offset1, index2, offset2, index3, offset_s
    integer :: interval
    real(dp) :: s, s2, w0(3), w1(3), w2(3)
    real(dp) :: dw0(3), dw1(3), dw2(3), wlist(100,3)
    integer :: k

    idx = find_segment(kernel, center, target)
    if (idx < 0) then
      print *, 'ERROR: segment not found center=', center, ' target=', target
      stop 1
    end if
    if (.not. kernel%segments(idx)%loaded) call load_segment_data(kernel, idx)

    init_e = kernel%segments(idx)%init_epoch
    intlen = kernel%segments(idx)%intlen
    n_int  = kernel%segments(idx)%n_intervals
    cc     = kernel%segments(idx)%coefficient_count

    call divmod_dp((tdb_whole - T0) * DAY_S - init_e, intlen, index1, offset1)
    call divmod_dp(tdb_frac * DAY_S, intlen, index2, offset2)
    call divmod_dp(offset1 + offset2, intlen, index3, offset_s)
    interval = int(index1 + index2 + index3)
    if (interval == n_int) then; interval = interval - 1; offset_s = offset_s + intlen; end if
    interval = interval + 1

    s = 2.0_dp * offset_s / intlen - 1.0_dp
    s2 = 2.0_dp * s

    w0 = 0.0_dp; w1 = 0.0_dp
    do k = 1, cc - 1
      w2 = w1; w1 = w0
      w0(1) = kernel%segments(idx)%coeffs(k,1,interval) + s2*w1(1) - w2(1)
      w0(2) = kernel%segments(idx)%coeffs(k,2,interval) + s2*w1(2) - w2(2)
      w0(3) = kernel%segments(idx)%coeffs(k,3,interval) + s2*w1(3) - w2(3)
      wlist(k,:) = w1
    end do
    pos(1) = kernel%segments(idx)%coeffs(cc,1,interval) + s*w0(1) - w1(1)
    pos(2) = kernel%segments(idx)%coeffs(cc,2,interval) + s*w0(2) - w1(2)
    pos(3) = kernel%segments(idx)%coeffs(cc,3,interval) + s*w0(3) - w1(3)

    dw0 = 0.0_dp; dw1 = 0.0_dp
    do k = 1, cc - 1
      dw2 = dw1; dw1 = dw0
      dw0(1) = 2.0_dp*wlist(k,1) + dw1(1)*s2 - dw2(1)
      dw0(2) = 2.0_dp*wlist(k,2) + dw1(2)*s2 - dw2(2)
      dw0(3) = 2.0_dp*wlist(k,3) + dw1(3)*s2 - dw2(3)
    end do
    vel(1) = w0(1) + s*dw0(1) - dw1(1)
    vel(2) = w0(2) + s*dw0(2) - dw1(2)
    vel(3) = w0(3) + s*dw0(3) - dw1(3)
    vel = vel / intlen * 2.0_dp * DAY_S
  end subroutine

  subroutine divmod_dp(a, b, q, r)
    real(dp), intent(in)  :: a, b
    real(dp), intent(out) :: q, r
    q = floor(a / b)
    r = a - q * b
  end subroutine

  subroutine spk_close(kernel)
    type(spk_kernel), intent(inout) :: kernel
    if (kernel%unit_num /= -1) close(kernel%unit_num)
    kernel%unit_num = -1
  end subroutine
end module spk_reader_mod

! ─── Module 3: Observation reader ────────────────────────────────────
module obs_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_OBS = 10000

  type :: observation
    real(dp) :: jd
    integer  :: body       ! 1=Sun, 2=Moon
    integer  :: event      ! 1=R, 2=S, 3=I
    real(dp) :: alt_obs    ! degrees
    real(dp) :: az_obs     ! degrees
    real(dp) :: alt_true   ! degrees
    real(dp) :: az_true    ! degrees
    real(dp) :: dist       ! AU
  end type

contains

  subroutine read_observations(filename, obs, n_obs, jd_max)
    character(len=*), intent(in)  :: filename
    type(observation), intent(out) :: obs(MAX_OBS)
    integer, intent(out) :: n_obs
    real(dp), intent(in) :: jd_max

    integer :: u, ios
    character(len=512) :: line
    character(1) :: body_c, event_c
    real(dp) :: jd, alt_o, az_o, alt_t, az_t, dist
    integer :: bi, ei

    open(newunit=u, file=filename, status='old', action='read')
    n_obs = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      read(line, *, iostat=ios) jd, body_c, event_c, alt_o, az_o, alt_t, az_t, dist
      if (ios /= 0) cycle
      if (jd > jd_max) cycle
      if (body_c == 'S') then; bi = 1; else; bi = 2; end if
      if (event_c == 'R') then; ei = 1
      else if (event_c == 'S') then; ei = 2
      else; ei = 3; end if
      n_obs = n_obs + 1
      if (n_obs > MAX_OBS) then
        n_obs = MAX_OBS; exit
      end if
      obs(n_obs)%jd       = jd
      obs(n_obs)%body     = bi
      obs(n_obs)%event    = ei
      obs(n_obs)%alt_obs  = alt_o
      obs(n_obs)%az_obs   = az_o
      obs(n_obs)%alt_true = alt_t
      obs(n_obs)%az_true  = az_t
      obs(n_obs)%dist     = dist
    end do
    close(u)
  end subroutine

end module obs_reader_mod

! ─── Module 4: N-body dynamics (Yoshida 4th-order) ──────────────────
module nbody2_mod
  use constants_mod
  implicit none
  integer, parameter :: NB = 2   ! Sun, Earth

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
    real(dp) :: r(3), r2, rinv3

    acc = 0.0_dp
    r = pos(2,:) - pos(1,:)
    r2 = dot_product(r, r)
    rinv3 = 1.0_dp / (r2 * sqrt(r2))
    acc(1,:) =  gm(2) * r * rinv3
    acc(2,:) = -gm(1) * r * rinv3
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

    ! 3-hour max sub-steps for stability
    dt_step = sign(min(abs(dt_total), 10800.0_dp), dt_total)
    n_steps = max(1, nint(abs(dt_total) / abs(dt_step)))
    dt_step = dt_total / real(n_steps, dp)

    do i = 1, n_steps
      call yoshida4_step_nb(pos, vel, gm, dt_step)
    end do
  end subroutine

end module nbody2_mod

! ─── Module 5: Keplerian elements + observation model ────────────────
module kepler_obs_mod
  use constants_mod
  use nbody2_mod, only: NB, propagate_nbody
  implicit none
  integer, parameter :: NS = 7   ! EKF state dimension
  real(dp), parameter :: GM_EARTH_KNOWN = 398600.435507_dp  ! fixed split

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
  subroutine predict_sun_altaz(x, t_epoch, jd_obs, lat_deg, lon_deg, &
                                alt_deg, az_deg)
    real(dp), intent(in)  :: x(NS), t_epoch, jd_obs, lat_deg, lon_deg
    real(dp), intent(out) :: alt_deg, az_deg

    real(dp) :: ecc, inc, raan_v, argp_v, M0, mu, a_val
    real(dp) :: dt_s
    real(dp) :: r_rel(3), v_rel(3), r_sun(3), v_earth(3)
    real(dp) :: pos_nb(NB,3), vel_nb(NB,3), gm_nb(NB), gm_s
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
    argp_v = x(4); M0 = x(5); mu = x(6); a_val = x(7)

    ! Convert elements at epoch to Cartesian (Earth relative to Sun)
    call kepler_to_cart(a_val, ecc, inc, raan_v, argp_v, M0, mu, r_rel, v_rel)

    ! Set up N-body initial conditions (barycentric frame)
    gm_s = mu - GM_EARTH_KNOWN
    pos_nb(1,:) = -(GM_EARTH_KNOWN / mu) * r_rel   ! Sun
    vel_nb(1,:) = -(GM_EARTH_KNOWN / mu) * v_rel
    pos_nb(2,:) =  (gm_s / mu) * r_rel              ! Earth
    vel_nb(2,:) =  (gm_s / mu) * v_rel
    gm_nb(1) = gm_s
    gm_nb(2) = GM_EARTH_KNOWN

    ! Propagate from epoch to observation time
    dt_s = (jd_obs - t_epoch) * DAY_S
    call propagate_nbody(pos_nb, vel_nb, gm_nb, dt_s)

    ! Sun geocentric position and Earth velocity
    r_sun   = pos_nb(1,:) - pos_nb(2,:)
    v_earth = vel_nb(2,:)

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
    v_obs = v_earth / C_LIGHT_KMS
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

end module kepler_obs_mod

! ═══════════════════════════════════════════════════════════════════════
!  Main program
! ═══════════════════════════════════════════════════════════════════════
program sun_ekf
  use constants_mod
  use spk_reader_mod
  use obs_reader_mod
  use kepler_obs_mod
  implicit none

  ! GM values (DE440)
  real(dp), parameter :: GM_SUN   = 1.32712440041279419d11
  real(dp), parameter :: GM_EARTH = 3.98600435507d5

  ! State vectors
  real(dp) :: x(NS), x_true(NS), x_init(NS)
  real(dp) :: t_epoch

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
  real(dp) :: corr_mat(NS,NS), sig_i, sig_j
  character(len=7) :: param_names(NS)

  ! ═══════════════════════════════════════════════════
  lat_obs = 40.0_dp
  lon_obs = 0.0_dp

  print '(A)', '════════════════════════════════════════════════════════'
  print '(A)', '  2-Body Sun-Earth EKF (N-body propagation)'
  print '(A,F8.4,A,F8.4)', '  Observer: lat=', lat_obs, ' lon=', lon_obs
  print '(A)', '════════════════════════════════════════════════════════'
  print '(A)', ''
  print '(A)', '  State vector: x = [e, i, Omega, omega, M0, mu, a]'
  print '(A)', '    e     = eccentricity (0 = circle, 1 = parabola)'
  print '(A)', '    i     = inclination of orbit to equator (deg)'
  print '(A)', '            (should recover obliquity ~23.44 deg)'
  print '(A)', '    Omega = longitude of ascending node (deg)'
  print '(A)', '            (where orbit crosses equator going north)'
  print '(A)', '    omega = argument of periapsis (deg)'
  print '(A)', '            (angle from ascending node to closest approach)'
  print '(A)', '    M0    = mean anomaly at epoch (deg)'
  print '(A)', '            (linearized orbital phase at reference time)'
  print '(A)', '    mu    = gravitational parameter GM_sun + GM_earth (km^3/s^2)'
  print '(A)', '    a     = semi-major axis of orbit (km)'

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

  ! True state vector
  x_true(1) = ecc_comp
  x_true(2) = inc_comp
  x_true(3) = raan_comp
  x_true(4) = argp_comp
  x_true(5) = M0_comp
  x_true(6) = mu_total
  x_true(7) = a_comp

  print '(/,A)',          '  True Keplerian elements (J2000 equatorial):'
  print '(A,F12.8)',      '    e     = ', x_true(1)
  print '(A,F10.5,A)',    '    i     = ', x_true(2) * RAD2DEG, ' deg (obliquity)'
  print '(A,F10.5,A)',    '    Omega = ', x_true(3) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    omega = ', x_true(4) * RAD2DEG, ' deg'
  print '(A,F10.5,A)',    '    M0    = ', x_true(5) * RAD2DEG, ' deg'
  print '(A,ES20.12,A)',  '    mu    = ', x_true(6), ' km^3/s^2'
  print '(A,F12.1,A)',    '    a     = ', x_true(7), ' km'

  ! ── 3. Perturb initial state ──
  x_init(1) = x_true(1) * 1.20_dp            ! e: +20%
  x_init(2) = x_true(2) + 2.0_dp * DEG2RAD   ! i: +2 deg
  x_init(3) = x_true(3) + 5.0_dp * DEG2RAD   ! Omega: +5 deg
  x_init(4) = x_true(4) + 5.0_dp * DEG2RAD   ! omega: +5 deg
  x_init(5) = x_true(5) + 3.0_dp * DEG2RAD   ! M0: +3 deg
  x_init(6) = x_true(6) * 1.005_dp           ! mu: +0.5%
  x_init(7) = x_true(7)                      ! a: correct value

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
  print '(A,F12.1,A)',          '    a     = ', x_init(7), ' km (correct)'

  ! ── 4. Initialize EKF ──
  x = x_init

  ! Initial covariance (1-sigma uncertainties on diagonal)
  P_cov = 0.0_dp
  P_cov(1,1) = (0.005_dp)**2             ! sigma_e = 0.005
  P_cov(2,2) = (3.0_dp * DEG2RAD)**2     ! sigma_i = 3 deg
  P_cov(3,3) = (8.0_dp * DEG2RAD)**2     ! sigma_Omega = 8 deg
  P_cov(4,4) = (8.0_dp * DEG2RAD)**2     ! sigma_omega = 8 deg
  P_cov(5,5) = (5.0_dp * DEG2RAD)**2     ! sigma_M0 = 5 deg
  P_cov(6,6) = (0.01_dp * mu_total)**2   ! sigma_mu = 1%
  P_cov(7,7) = (0.01_dp * a_comp)**2     ! sigma_a = 1% (~1.5M km)

  print '(/,A)', '  Initial 1-sigma uncertainties:'
  print '(A,F10.6)',     '    sigma_e     = ', sqrt(P_cov(1,1))
  print '(A,F8.4,A)',    '    sigma_i     = ', sqrt(P_cov(2,2))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',    '    sigma_Omega = ', sqrt(P_cov(3,3))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',    '    sigma_omega = ', sqrt(P_cov(4,4))*RAD2DEG, ' deg'
  print '(A,F8.4,A)',    '    sigma_M0    = ', sqrt(P_cov(5,5))*RAD2DEG, ' deg'
  print '(A,ES10.3,A,F6.3,A)', '    sigma_mu    = ', sqrt(P_cov(6,6)), &
       '  (', sqrt(P_cov(6,6))/mu_total*100.0_dp, '%)'
  print '(A,ES10.3,A,F6.3,A)', '    sigma_a     = ', sqrt(P_cov(7,7)), &
       '  (', sqrt(P_cov(7,7))/a_comp*100.0_dp, '%)'

  Q_noise = 0.0_dp   ! elements are constant in 2-body

  sigma_obs = 60.0_dp / 3600.0_dp   ! 60 arcsec in degrees
  R_noise = 0.0_dp
  R_noise(1,1) = sigma_obs**2
  R_noise(2,2) = sigma_obs**2

  print '(A,F6.1,A)',    '    sigma_obs   = ', sigma_obs * 3600.0_dp, &
       ' arcsec (measurement noise, alt & az)'

  ! FD perturbation sizes
  delta(1) = 1.0d-7             ! e
  delta(2) = 1.0d-7             ! i (rad)
  delta(3) = 1.0d-7             ! Omega (rad)
  delta(4) = 1.0d-7             ! omega (rad)
  delta(5) = 1.0d-7             ! M0 (rad)
  delta(6) = mu_total * 1.0d-7  ! mu
  delta(7) = a_comp * 1.0d-7    ! a

  I_mat = 0.0_dp
  do i = 1, NS
    I_mat(i,i) = 1.0_dp
  end do

  ! ── 5. EKF loop ──
  print '(/,A)', '  Running EKF...'
  print '(A)', ''
  print '(A)', '  Column definitions:'
  print '(A)', '    Obs#     = observation number (cumulative)'
  print '(A)', '    cumRMSa" = cumulative RMS of altitude residuals (arcsec)'
  print '(A)', '    cumRMSz" = cumulative RMS of azimuth residuals (arcsec)'
  print '(A)', '    winRMSa" = windowed RMS of altitude residuals, last 200 obs (arcsec)'
  print '(A)', '    winRMSz" = windowed RMS of azimuth residuals, last 200 obs (arcsec)'
  print '(A)', '    e_err%   = relative error in eccentricity vs true (%)'
  print '(A)', '    i_err(d) = error in inclination vs true (degrees)'
  print '(A)', '    mu_err%  = relative error in mu vs true (%)'
  print '(A)', '    a_err%   = relative error in semi-major axis vs true (%)'
  print '(A)', '  ──────────────────────────────────────────────────────────────────────────'
  print '(A)', '   Obs#  cumRMSa" cumRMSz"  winRMSa"  winRMSz"  e_err%   i_err(d) mu_err%  a_err%'

  rms_alt = 0.0_dp
  rms_az  = 0.0_dp
  win_alt = 0.0_dp
  win_az  = 0.0_dp
  n_proc  = 0

  do k = 1, n_sun
    ! ── Predict (trivial: state unchanged, covariance grows by Q) ──
    P_pred = P_cov + Q_noise

    ! ── Predicted observation ──
    call predict_sun_altaz(x, t_epoch, s_jd(k), lat_obs, lon_obs, &
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
      call predict_sun_altaz(x_pert, t_epoch, s_jd(k), &
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
    if (x(7) < 0.0_dp) x(7) = a_comp * 0.9_dp

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
      print '(I7,2F10.1,2F10.1,F9.4,F10.5,F8.3,F8.3)', &
           n_proc, &
           sqrt(rms_alt / n_proc) * 3600.0_dp, &
           sqrt(rms_az  / n_proc) * 3600.0_dp, &
           sqrt(win_alt / j) * 3600.0_dp, &
           sqrt(win_az  / j) * 3600.0_dp, &
           (x(1)/x_true(1) - 1.0_dp) * 100.0_dp, &
           (x(2) - x_true(2)) * RAD2DEG, &
           (x(6)/x_true(6) - 1.0_dp) * 100.0_dp, &
           (x(7)/x_true(7) - 1.0_dp) * 100.0_dp
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
  print '(A,F12.1,A,F12.1,A,F8.4,A)', &
       '  a       ', x(7), ' km  true: ', x_true(7), ' km  err: ', &
       (x(7)/x_true(7)-1.0_dp)*100.0_dp, '%'

  print '(/,A)', '  Derived quantities:'
  print '(A,F10.5,A,F10.5,A)', '    Obliquity  = ', x(2)*RAD2DEG, &
       ' deg  (true: ', x_true(2)*RAD2DEG, ' deg)'
  print '(A,F12.3,A)', '    Period     = ', &
       TAU / sqrt(x(6) / x(7)**3) / DAY_S, ' days'
  print '(A,F12.3,A)', '    (true)     = ', &
       TAU / sqrt(x_true(6) / x_true(7)**3) / DAY_S, ' days'

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
  print '(A,ES10.3,A,F8.5,A)', '    sigma_a     = ', &
       sqrt(max(0.0_dp, P_cov(7,7))), &
       '  (', sqrt(max(0.0_dp, P_cov(7,7)))/x(7)*100.0_dp, '%)'

  ! Correlation matrix
  param_names(1) = '  e    '
  param_names(2) = '  i    '
  param_names(3) = '  Omega'
  param_names(4) = '  omega'
  param_names(5) = '  M0   '
  param_names(6) = '  mu   '
  param_names(7) = '  a    '
  do i = 1, NS
    sig_i = sqrt(max(0.0_dp, P_cov(i,i)))
    do j = 1, NS
      sig_j = sqrt(max(0.0_dp, P_cov(j,j)))
      if (sig_i > 0.0_dp .and. sig_j > 0.0_dp) then
        corr_mat(i,j) = P_cov(i,j) / (sig_i * sig_j)
      else
        corr_mat(i,j) = 0.0_dp
      end if
    end do
  end do

  print '(/,A)', '  Correlation matrix:'
  print '(A,7(A7,1X))', '          ', (param_names(j), j=1,NS)
  do i = 1, NS
    print '(A,7F8.4)', param_names(i), (corr_mat(i,j), j=1,NS)
  end do

  print '(/,A)', '  Note: "True" elements are osculating at epoch. The EKF'
  print '(A)',   '  estimates best-fit MEAN elements over the full arc.'
  print '(A)',   '  The e difference (~5%) is expected: planetary perturbations'
  print '(A)',   '  make the osculating e vary over the year.'
  print '(A)', '════════════════════════════════════════════════════════'

end program sun_ekf
