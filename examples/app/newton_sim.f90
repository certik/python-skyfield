module nbody_mod
  use constants_mod
  implicit none

  integer, parameter :: NBODIES = 10

  real(dp), parameter :: GM(NBODIES) = [ &
    132712440041.279419_dp, &
    22031.868551_dp, &
    324858.592000_dp, &
    398600.435507_dp, &
    4902.800118_dp, &
    42828.375816_dp, &
    126712764.100000_dp, &
    37940584.841800_dp, &
    5794556.400000_dp, &
    6836527.100580_dp ]

  real(dp), parameter :: W0 = -(2.0_dp ** (1.0_dp / 3.0_dp)) / &
      (2.0_dp - 2.0_dp ** (1.0_dp / 3.0_dp))
  real(dp), parameter :: W1 = 1.0_dp / (2.0_dp - 2.0_dp ** (1.0_dp / 3.0_dp))
  real(dp), parameter :: C1 = 0.5_dp * W1
  real(dp), parameter :: C4 = C1
  real(dp), parameter :: C2 = 0.5_dp * (W0 + W1)
  real(dp), parameter :: C3 = C2
  real(dp), parameter :: D1 = W1
  real(dp), parameter :: D2 = W0
  real(dp), parameter :: D3 = W1

contains

  subroutine compute_accelerations(pos, acc)
    real(dp), intent(in)  :: pos(NBODIES,3)
    real(dp), intent(out) :: acc(NBODIES,3)
    integer :: i, j
    real(dp) :: r(3), r2, rinv3

    acc = 0.0_dp
    do i = 1, NBODIES - 1
      do j = i + 1, NBODIES
        r = pos(j,:) - pos(i,:)
        r2 = dot_product(r, r)
        rinv3 = 1.0_dp / (r2 * sqrt(r2))
        acc(i,:) = acc(i,:) + GM(j) * r * rinv3
        acc(j,:) = acc(j,:) - GM(i) * r * rinv3
      end do
    end do
  end subroutine

  subroutine yoshida4_step(pos, vel, dt)
    real(dp), intent(inout) :: pos(NBODIES,3), vel(NBODIES,3)
    real(dp), intent(in)    :: dt
    real(dp) :: acc(NBODIES,3)

    pos = pos + C1 * dt * vel
    call compute_accelerations(pos, acc)
    vel = vel + D1 * dt * acc

    pos = pos + C2 * dt * vel
    call compute_accelerations(pos, acc)
    vel = vel + D2 * dt * acc

    pos = pos + C3 * dt * vel
    call compute_accelerations(pos, acc)
    vel = vel + D3 * dt * acc

    pos = pos + C4 * dt * vel
  end subroutine
end module nbody_mod


module chebyshev_fit_mod
  use constants_mod
  implicit none
contains

  subroutine fit_chebyshev(x, y, ncoeff, coeffs)
    real(dp), intent(in) :: x(0:), y(0:)
    integer, intent(in)  :: ncoeff
    real(dp), intent(out) :: coeffs(ncoeff)
    integer :: n_samples, i, j, k
    real(dp), allocatable :: normal(:,:), rhs(:), basis(:), sol(:)

    n_samples = ubound(x, 1)
    if (ubound(y, 1) /= n_samples) then
      print *, 'ERROR: fit_chebyshev input lengths mismatch'
      stop 1
    end if

    allocate(normal(ncoeff,ncoeff), rhs(ncoeff), basis(ncoeff), sol(ncoeff))
    normal = 0.0_dp
    rhs = 0.0_dp

    do k = 0, n_samples
      call chebyshev_basis(x(k), ncoeff, basis)
      do i = 1, ncoeff
        rhs(i) = rhs(i) + basis(i) * y(k)
        do j = 1, ncoeff
          normal(i,j) = normal(i,j) + basis(i) * basis(j)
        end do
      end do
    end do

    call solve_linear_system(normal, rhs, ncoeff, sol)
    coeffs = sol
    deallocate(normal, rhs, basis, sol)
  end subroutine

  subroutine chebyshev_basis(x, ncoeff, basis)
    real(dp), intent(in) :: x
    integer, intent(in) :: ncoeff
    real(dp), intent(out) :: basis(ncoeff)
    integer :: i

    basis(1) = 1.0_dp
    if (ncoeff >= 2) basis(2) = x
    do i = 3, ncoeff
      basis(i) = 2.0_dp * x * basis(i-1) - basis(i-2)
    end do
  end subroutine

  subroutine solve_linear_system(a, b, n, x)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n,n)
    real(dp), intent(inout) :: b(n)
    real(dp), intent(out) :: x(n)
    integer :: i, j, k, pivot
    real(dp) :: maxabs, factor, temp_b
    real(dp) :: temp_row(n)

    do k = 1, n - 1
      pivot = k
      maxabs = abs(a(k,k))
      do i = k + 1, n
        if (abs(a(i,k)) > maxabs) then
          maxabs = abs(a(i,k))
          pivot = i
        end if
      end do
      if (maxabs <= tiny(1.0_dp)) then
        print *, 'ERROR: singular normal matrix in Chebyshev fit'
        stop 1
      end if

      if (pivot /= k) then
        temp_row = a(k,:)
        a(k,:) = a(pivot,:)
        a(pivot,:) = temp_row
        temp_b = b(k)
        b(k) = b(pivot)
        b(pivot) = temp_b
      end if

      do i = k + 1, n
        factor = a(i,k) / a(k,k)
        a(i,k:n) = a(i,k:n) - factor * a(k,k:n)
        b(i) = b(i) - factor * b(k)
      end do
    end do

    if (abs(a(n,n)) <= tiny(1.0_dp)) then
      print *, 'ERROR: singular matrix at back substitution'
      stop 1
    end if

    x(n) = b(n) / a(n,n)
    do i = n - 1, 1, -1
      x(i) = (b(i) - sum(a(i,i+1:n) * x(i+1:n))) / a(i,i)
    end do
  end subroutine
end module chebyshev_fit_mod


module spk_writer_mod
  use constants_mod
  implicit none

  type :: spk_out_segment
    integer :: target
    integer :: center
    integer :: frame = 1
    integer :: data_type = 2
    integer :: ncoeff
    integer :: n_intervals
    real(dp) :: start_second
    real(dp) :: end_second
    real(dp) :: intlen
    real(dp), allocatable :: coeff(:,:,:)  ! (ncoeff, 3, n_intervals), low→high
  end type

contains

  subroutine spk_create(filename, segments, nseg)
    character(len=*), intent(in) :: filename
    type(spk_out_segment), intent(in) :: segments(nseg)
    integer, intent(in) :: nseg
    integer :: start_i(nseg), end_i(nseg)
    integer :: word_ptr, last_data_word
    integer :: summary_rec, name_rec, free_word
    integer :: i, rec, k, rsize
    integer :: base_pos, summary_pos, name_pos
    integer :: u
    real(dp) :: mid, radius
    character(8) :: locidw, locfmt
    character(60) :: locifn
    character(603) :: prenul
    character(28) :: ftpstr
    character(297) :: pstnul
    character(40) :: seg_name
    character(1024) :: blank_record

    word_ptr = 257
    do i = 1, nseg
      rsize = 2 + 3 * segments(i)%ncoeff
      start_i(i) = word_ptr
      end_i(i) = word_ptr + rsize * segments(i)%n_intervals + 3
      word_ptr = end_i(i) + 1
    end do

    last_data_word = word_ptr - 1
    summary_rec = (last_data_word - 1) / 128 + 2
    name_rec = summary_rec + 1
    free_word = name_rec * 128 + 1

    open(newunit=u, file=filename, access='stream', form='unformatted', &
         status='replace', action='write')

    blank_record = repeat(' ', 1024)
    write(u, pos=1) blank_record
    write(u, pos=1025) blank_record
    write(u, pos=(summary_rec - 1) * 1024 + 1) blank_record
    write(u, pos=(name_rec - 1) * 1024 + 1) blank_record

    locidw = 'NAIF/DAF'
    locfmt = 'LTL-IEEE'
    locifn = repeat(' ', 60)
    locifn(1:min(60, len_trim(filename))) = filename(1:min(60, len_trim(filename)))
    prenul = repeat(' ', 603)
    ftpstr = repeat(' ', 28)
    ftpstr(1:18) = 'FTPSTR:TEST:ENDFTP'
    pstnul = repeat(' ', 297)

    write(u, pos=1) locidw
    write(u, pos=9) 2
    write(u, pos=13) 6
    write(u, pos=17) locifn
    write(u, pos=77) summary_rec
    write(u, pos=81) summary_rec
    write(u, pos=85) free_word
    write(u, pos=89) locfmt
    write(u, pos=97) prenul
    write(u, pos=700) ftpstr
    write(u, pos=728) pstnul

    do i = 1, nseg
      rsize = 2 + 3 * segments(i)%ncoeff
      do rec = 1, segments(i)%n_intervals
        mid = segments(i)%start_second + (real(rec, dp) - 0.5_dp) * segments(i)%intlen
        radius = 0.5_dp * segments(i)%intlen

        base_pos = (start_i(i) - 1 + (rec - 1) * rsize) * 8 + 1
        write(u, pos=base_pos) mid, radius
        write(u, pos=base_pos + 16) &
            (segments(i)%coeff(k,1,rec), k = 1, segments(i)%ncoeff), &
            (segments(i)%coeff(k,2,rec), k = 1, segments(i)%ncoeff), &
            (segments(i)%coeff(k,3,rec), k = 1, segments(i)%ncoeff)
      end do

      base_pos = (end_i(i) - 4) * 8 + 1
      write(u, pos=base_pos) segments(i)%start_second, segments(i)%intlen, &
          real(rsize, dp), real(segments(i)%n_intervals, dp)
    end do

    base_pos = (summary_rec - 1) * 1024 + 1
    write(u, pos=base_pos) 0.0_dp, 0.0_dp, real(nseg, dp)
    do i = 1, nseg
      summary_pos = base_pos + 24 + (i - 1) * 40
      write(u, pos=summary_pos) segments(i)%start_second, segments(i)%end_second, &
          segments(i)%target, segments(i)%center, segments(i)%frame, &
          segments(i)%data_type, start_i(i), end_i(i)
    end do

    base_pos = (name_rec - 1) * 1024 + 1
    do i = 1, nseg
      seg_name = repeat(' ', 40)
      write(seg_name, '(A,I0,A,I0,A)') 'NEWTON ', segments(i)%center, '->', &
          segments(i)%target, ' TYPE2'
      name_pos = base_pos + (i - 1) * 40
      write(u, pos=name_pos) seg_name
    end do

    close(u)
  end subroutine
end module spk_writer_mod


program newton_sim
  use constants_mod
  use spk_reader_mod
  use nbody_mod
  use chebyshev_fit_mod
  use spk_writer_mod
  implicit none

  integer, parameter :: NSEG = 11
  integer, parameter :: TARGETS(NSEG) = [1,2,3,4,5,6,7,8,10,301,399]
  integer, parameter :: CENTERS(NSEG) = [0,0,0,0,0,0,0,0,0,3,3]
  integer, parameter :: INTERVAL_DAYS(NSEG) = [8,16,16,32,32,32,32,32,16,4,4]
  integer, parameter :: NCOEFFS(NSEG) = [14,10,13,11,8,7,6,6,11,13,13]

  real(dp), parameter :: START_JD = 2415020.5_dp
  real(dp), parameter :: TARGET_END_JD = 2462502.5_dp
  real(dp), parameter :: DT_DAY = 0.125_dp

  type(spk_kernel) :: source_kernel, ref_kernel, newton_kernel
  type(spk_out_segment), allocatable :: segments(:)

  integer :: n_steps, step, s, interval, comp, k
  integer :: n_intervals(NSEG), n_eval, max_coeff, max_eval
  real(dp) :: start_second, target_span_days, dt_sec
  real(dp) :: interval_seconds(NSEG), segment_end_second(NSEG), global_end_second
  real(dp), allocatable :: traj(:,:,:), seg_traj(:,:,:)
  real(dp) :: pos(NBODIES,3), vel(NBODIES,3)
  real(dp) :: earth(3), moon(3), emb(3)
  real(dp), allocatable :: x_nodes(:), y_nodes(:), coeff(:)
  real(dp) :: node_second
  real(dp) :: p_ref(3), p_new(3), emb_ref(3), emb_new(3), mrel_ref(3), mrel_new(3)
  real(dp) :: sun_err_km, moon_err_km

  target_span_days = TARGET_END_JD - START_JD
  dt_sec = DT_DAY * DAY_S
  start_second = (START_JD - T0) * DAY_S

  do s = 1, NSEG
    interval_seconds(s) = real(INTERVAL_DAYS(s), dp) * DAY_S
    n_intervals(s) = ceiling(target_span_days / real(INTERVAL_DAYS(s), dp))
    segment_end_second(s) = start_second + real(n_intervals(s), dp) * interval_seconds(s)
  end do

  global_end_second = maxval(segment_end_second)
  n_steps = nint((global_end_second - start_second) / dt_sec)

  print '(A,F10.1,A)', 'Integrating span: ', (global_end_second - start_second) / DAY_S, ' days'
  print '(A,I0)', 'Integrator steps: ', n_steps

  allocate(traj(NBODIES,3,0:n_steps))
  allocate(seg_traj(NSEG,3,0:n_steps))

  call spk_open('de440s.bsp', source_kernel)
  call load_initial_conditions(source_kernel, START_JD, pos, vel)
  call spk_close(source_kernel)

  traj(:,:,0) = pos
  do step = 1, n_steps
    call yoshida4_step(pos, vel, dt_sec)
    traj(:,:,step) = pos
  end do

  do step = 0, n_steps
    earth = traj(4,:,step)
    moon = traj(5,:,step)
    emb = (GM(4) * earth + GM(5) * moon) / (GM(4) + GM(5))

    seg_traj(1,:,step) = traj(2,:,step)     ! Mercury wrt SSB
    seg_traj(2,:,step) = traj(3,:,step)     ! Venus wrt SSB
    seg_traj(3,:,step) = emb                ! EMB wrt SSB
    seg_traj(4,:,step) = traj(6,:,step)     ! Mars wrt SSB
    seg_traj(5,:,step) = traj(7,:,step)     ! Jupiter wrt SSB
    seg_traj(6,:,step) = traj(8,:,step)     ! Saturn wrt SSB
    seg_traj(7,:,step) = traj(9,:,step)     ! Uranus wrt SSB
    seg_traj(8,:,step) = traj(10,:,step)    ! Neptune wrt SSB
    seg_traj(9,:,step) = traj(1,:,step)     ! Sun wrt SSB
    seg_traj(10,:,step) = moon - emb        ! Moon wrt EMB
    seg_traj(11,:,step) = earth - emb       ! Earth wrt EMB
  end do

  allocate(segments(NSEG))
  do s = 1, NSEG
    segments(s)%target = TARGETS(s)
    segments(s)%center = CENTERS(s)
    segments(s)%frame = 1
    segments(s)%data_type = 2
    segments(s)%ncoeff = NCOEFFS(s)
    segments(s)%n_intervals = n_intervals(s)
    segments(s)%start_second = start_second
    segments(s)%end_second = segment_end_second(s)
    segments(s)%intlen = interval_seconds(s)
    allocate(segments(s)%coeff(NCOEFFS(s),3,n_intervals(s)))
  end do

  max_coeff = maxval(NCOEFFS)
  max_eval = 2 * max_coeff
  allocate(x_nodes(0:max_eval), y_nodes(0:max_eval), coeff(max_coeff))

  do s = 1, NSEG
    n_eval = 2 * NCOEFFS(s)
    do interval = 1, n_intervals(s)
      do comp = 1, 3
        do k = 0, n_eval
          x_nodes(k) = cos(PI * real(k, dp) / real(n_eval, dp))
          node_second = start_second + (real(interval, dp) - 0.5_dp) * interval_seconds(s) + &
              0.5_dp * interval_seconds(s) * x_nodes(k)
          y_nodes(k) = sample_series(seg_traj(s,comp,:), node_second - start_second, dt_sec)
        end do
        call fit_chebyshev(x_nodes(0:n_eval), y_nodes(0:n_eval), NCOEFFS(s), coeff(1:NCOEFFS(s)))
        segments(s)%coeff(:,comp,interval) = coeff(1:NCOEFFS(s))
      end do
    end do
  end do

  call spk_create('newton.bsp', segments, NSEG)
  print '(A)', 'Wrote newton.bsp'

  call spk_open('de440s.bsp', ref_kernel)
  call spk_open('newton.bsp', newton_kernel)

  call spk_compute(ref_kernel, 0, 10, 2460676.5_dp, 0.0_dp, p_ref)
  call spk_compute(newton_kernel, 0, 10, 2460676.5_dp, 0.0_dp, p_new)
  sun_err_km = norm3(p_new - p_ref)

  call spk_compute(ref_kernel, 0, 3, 2460676.5_dp, 0.0_dp, emb_ref)
  call spk_compute(ref_kernel, 3, 301, 2460676.5_dp, 0.0_dp, mrel_ref)
  call spk_compute(newton_kernel, 0, 3, 2460676.5_dp, 0.0_dp, emb_new)
  call spk_compute(newton_kernel, 3, 301, 2460676.5_dp, 0.0_dp, mrel_new)
  moon_err_km = norm3((emb_new + mrel_new) - (emb_ref + mrel_ref))

  print '(A,ES12.4,A)', 'Sun position error at 2025-01-01:  ', sun_err_km, ' km'
  print '(A,ES12.4,A)', 'Moon position error at 2025-01-01: ', moon_err_km, ' km'

  call spk_close(ref_kernel)
  call spk_close(newton_kernel)

contains

  subroutine load_initial_conditions(kernel, jd, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in) :: jd
    real(dp), intent(out) :: pos(NBODIES,3), vel(NBODIES,3)
    real(dp) :: p1(3), p2(3), v1(3), v2(3)

    call read_state(kernel, jd, 0, 10, pos(1,:), vel(1,:))

    call read_state(kernel, jd, 0, 1, p1, v1)
    call read_state(kernel, jd, 1, 199, p2, v2)
    pos(2,:) = p1 + p2
    vel(2,:) = v1 + v2

    call read_state(kernel, jd, 0, 2, p1, v1)
    call read_state(kernel, jd, 2, 299, p2, v2)
    pos(3,:) = p1 + p2
    vel(3,:) = v1 + v2

    call read_state(kernel, jd, 0, 3, p1, v1)
    call read_state(kernel, jd, 3, 399, p2, v2)
    pos(4,:) = p1 + p2
    vel(4,:) = v1 + v2

    call read_state(kernel, jd, 0, 3, p1, v1)
    call read_state(kernel, jd, 3, 301, p2, v2)
    pos(5,:) = p1 + p2
    vel(5,:) = v1 + v2

    call read_state(kernel, jd, 0, 4, pos(6,:), vel(6,:))
    call read_state(kernel, jd, 0, 5, pos(7,:), vel(7,:))
    call read_state(kernel, jd, 0, 6, pos(8,:), vel(8,:))
    call read_state(kernel, jd, 0, 7, pos(9,:), vel(9,:))
    call read_state(kernel, jd, 0, 8, pos(10,:), vel(10,:))
  end subroutine

  subroutine read_state(kernel, jd, center, target, pos, vel)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in) :: jd
    integer, intent(in) :: center, target
    real(dp), intent(out) :: pos(3), vel(3)
    real(dp) :: vel_day(3)
    call spk_compute_and_diff(kernel, center, target, jd, 0.0_dp, pos, vel_day)
    vel = vel_day / DAY_S
  end subroutine

  function sample_series(series, t_rel_second, dt_second) result(value)
    real(dp), intent(in) :: series(0:)
    real(dp), intent(in) :: t_rel_second, dt_second
    real(dp) :: value
    integer :: nmax, i0
    real(dp) :: idx, frac

    nmax = ubound(series, 1)
    if (t_rel_second <= 0.0_dp) then
      value = series(0)
      return
    end if
    if (t_rel_second >= real(nmax, dp) * dt_second) then
      value = series(nmax)
      return
    end if

    idx = t_rel_second / dt_second
    i0 = int(floor(idx))
    if (i0 >= nmax) then
      value = series(nmax)
      return
    end if
    frac = idx - real(i0, dp)
    value = (1.0_dp - frac) * series(i0) + frac * series(i0 + 1)
  end function

  function norm3(v) result(r)
    real(dp), intent(in) :: v(3)
    real(dp) :: r
    r = sqrt(dot_product(v, v))
  end function
end program newton_sim
