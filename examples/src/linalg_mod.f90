! ─────────────────────────────────────────────────────────────────────
!  linalg_mod — small linear-algebra helpers (3×3 matrices, 3-vectors)
! ─────────────────────────────────────────────────────────────────────
module linalg_mod
  use constants_mod
  implicit none
contains

  ! ── 3×3 matrix multiply ──
  function mat33_mul(A, B) result(C)
    real(dp), intent(in) :: A(3,3), B(3,3)
    real(dp) :: C(3,3)
    integer :: i, j, k
    C = 0.0_dp
    do j = 1, 3
      do k = 1, 3
        do i = 1, 3
          C(i,j) = C(i,j) + A(i,k) * B(k,j)
        end do
      end do
    end do
  end function

  ! ── 3×3 × 3-vector ──
  function mat33_vec(A, x) result(y)
    real(dp), intent(in) :: A(3,3), x(3)
    real(dp) :: y(3)
    integer :: i, k
    y = 0.0_dp
    do k = 1, 3
      do i = 1, 3
        y(i) = y(i) + A(i,k) * x(k)
      end do
    end do
  end function

  ! ── transpose ──
  function mat33_T(A) result(AT)
    real(dp), intent(in) :: A(3,3)
    real(dp) :: AT(3,3)
    integer :: i, j
    do i = 1, 3
      do j = 1, 3
        AT(i,j) = A(j,i)
      end do
    end do
  end function

  ! ── 3-vector dot ──
  function dot3(a, b) result(d)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: d
    d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function

  ! ── 3-vector length ──
  function vec_len(v) result(r)
    real(dp), intent(in) :: v(3)
    real(dp) :: r
    r = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
  end function

  ! ── rotation matrices ──
  function rot_x(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [1.0_dp,  0.0_dp, 0.0_dp]
    R(2,:) = [0.0_dp,  c,      -s     ]
    R(3,:) = [0.0_dp,  s,       c     ]
  end function

  function rot_y(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [ c,      0.0_dp, s     ]
    R(2,:) = [ 0.0_dp, 1.0_dp, 0.0_dp]
    R(3,:) = [-s,      0.0_dp, c     ]
  end function

  function rot_z(angle) result(R)
    real(dp), intent(in) :: angle
    real(dp) :: R(3,3), c, s
    c = cos(angle); s = sin(angle)
    R(1,:) = [ c, -s,      0.0_dp]
    R(2,:) = [ s,  c,      0.0_dp]
    R(3,:) = [ 0.0_dp, 0.0_dp, 1.0_dp]
  end function

  ! ── to_spherical: xyz → (r, elevation, azimuth) ──
  subroutine to_spherical(xyz, r, elev, azim)
    real(dp), intent(in)  :: xyz(3)
    real(dp), intent(out) :: r, elev, azim
    real(dp) :: eps
    eps = tiny(1.0_dp)
    r = vec_len(xyz)
    elev = asin(xyz(3) / (r + eps))
    azim = mod(atan2(xyz(2), xyz(1)), TAU)
    if (azim < 0.0_dp) azim = azim + TAU
  end subroutine
end module linalg_mod
