! ─────────────────────────────────────────────────────────────────────
!  nbody2_mod — 2-body Yoshida 4th-order symplectic integrator
! ─────────────────────────────────────────────────────────────────────
module nbody2_mod
  use constants_mod
  implicit none
  integer, parameter :: NB = 2   ! Sun, Earth (or Earth, Moon)

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
