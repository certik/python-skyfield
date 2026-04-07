! ─────────────────────────────────────────────────────────────────────
!  obs_reader_mod — reads synthetic Sun/Moon alt/az observations
! ─────────────────────────────────────────────────────────────────────
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
