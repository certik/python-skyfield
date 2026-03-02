! ═══════════════════════════════════════════════════════════════════════
!  newton_sim.f90
!
!  Build a Newtonian SPK kernel (newton.bsp) by:
!    1) Reading DE440s initial states at J2010.0
!    2) Integrating a 10-body Newtonian model with Yoshida 4th-order scheme
!    3) Fitting type-2 Chebyshev segments
!    4) Writing a minimal DAF/SPK file
! ═══════════════════════════════════════════════════════════════════════

module constants_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI    = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: T0    = 2451545.0_dp
  real(dp), parameter :: DAY_S = 86400.0_dp
end module constants_mod


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
    real(dp) :: init_epoch
    real(dp) :: intlen
    integer  :: n_intervals
    integer  :: coefficient_count
    integer  :: component_count
    real(dp), allocatable :: coeffs(:,:,:)  ! (coeff, comp, interval), reversed
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
        if (seg_idx > MAX_SEGMENTS) then
          print *, 'ERROR: too many segments in SPK'
          stop 1
        end if

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
        idx = i
        return
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
    read(u, pos=pos) meta(1), meta(2), meta(3), meta(4)

    kernel%segments(idx)%init_epoch = meta(1)
    kernel%segments(idx)%intlen     = meta(2)
    rsize_i = int(meta(3))
    n_i     = int(meta(4))
    kernel%segments(idx)%n_intervals = n_i

    if (seg%data_type == 2) then
      comp_count = 3
    else
      comp_count = 6
    end if
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
              raw((rec - 1) * rsize_i + 2 + (c - 1) * coeff_count + k)
        end do
      end do
    end do

    deallocate(raw)
    kernel%segments(idx)%loaded = .true.
  end subroutine

  subroutine spk_compute(kernel, center, target, tdb_whole, tdb_frac, pos)
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in) :: center, target
    real(dp), intent(in) :: tdb_whole, tdb_frac
    real(dp), intent(out) :: pos(3)
    real(dp) :: vel(3)
    call spk_compute_and_diff(kernel, center, target, tdb_whole, tdb_frac, pos, vel)
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
    real(dp) :: s, s2
    real(dp) :: w0(3), w1(3), w2(3)
    real(dp) :: dw0(3), dw1(3), dw2(3)
    real(dp) :: wlist(100,3)
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

    if (interval == n_int) then
      interval = interval - 1
      offset_s = offset_s + intlen
    end if
    interval = interval + 1

    s = 2.0_dp * offset_s / intlen - 1.0_dp
    s2 = 2.0_dp * s

    w0 = 0.0_dp
    w1 = 0.0_dp
    do k = 1, cc - 1
      w2 = w1
      w1 = w0
      w0(1) = kernel%segments(idx)%coeffs(k,1,interval) + s2 * w1(1) - w2(1)
      w0(2) = kernel%segments(idx)%coeffs(k,2,interval) + s2 * w1(2) - w2(2)
      w0(3) = kernel%segments(idx)%coeffs(k,3,interval) + s2 * w1(3) - w2(3)
      wlist(k,:) = w1
    end do

    pos(1) = kernel%segments(idx)%coeffs(cc,1,interval) + s * w0(1) - w1(1)
    pos(2) = kernel%segments(idx)%coeffs(cc,2,interval) + s * w0(2) - w1(2)
    pos(3) = kernel%segments(idx)%coeffs(cc,3,interval) + s * w0(3) - w1(3)

    dw0 = 0.0_dp
    dw1 = 0.0_dp
    do k = 1, cc - 1
      dw2 = dw1
      dw1 = dw0
      dw0(1) = 2.0_dp * wlist(k,1) + dw1(1) * s2 - dw2(1)
      dw0(2) = 2.0_dp * wlist(k,2) + dw1(2) * s2 - dw2(2)
      dw0(3) = 2.0_dp * wlist(k,3) + dw1(3) * s2 - dw2(3)
    end do

    vel(1) = w0(1) + s * dw0(1) - dw1(1)
    vel(2) = w0(2) + s * dw0(2) - dw1(2)
    vel(3) = w0(3) + s * dw0(3) - dw1(3)
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
