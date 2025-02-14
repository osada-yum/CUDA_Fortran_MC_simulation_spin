!> CUDAFortran implementation for the 2-dimensional XY model with periodic boundary condition.
module xy2d_periodic_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), public, protected :: xy2d_gpu_stat

  integer(int64), parameter :: NUM_THREADS = 512

  public :: xy2d_gpu
  type :: xy2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     real(real64), allocatable, device :: spins_(:, :, :)
     real(real64), allocatable, device :: autocorrelation_start_spins_(:, :, :)
     real(real64), allocatable, device :: randoms_(:, :)
     real(real64), allocatable, device :: candidates_(:, :)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_xy2d_gpu
     procedure, pass :: skip_curand => skip_curand_xy2d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_xy2d_gpu
     procedure, pass :: set_random_spin => set_random_spin_xy2d_gpu
     procedure, pass :: set_random_small_spin => set_random_small_spin_xy2d_gpu
     procedure, pass :: set_random_near_spin => set_random_near_spin_xy2d_gpu
     procedure, pass :: set_finite_magne_spin => set_finite_magne_spin_xy2d_gpu
     procedure, pass :: rotate_summation_magne_toward_xaxis => rotate_summation_magne_toward_xaxis_xy2d_gpu
     procedure, pass :: rotate_summation_magne_and_autocorrelation_toward_xaxis => &
          & rotate_summation_magne_and_autocorrelation_toward_xaxis_xy2d_gpu
     ! procedure, pass :: rotate_summation_magne_and_autocorrelation_toward_xaxis_updown_randomly => &
     !      & rotate_summation_magne_and_autocorrelation_toward_xaxis_updown_randomly_xy2d_gpu
     procedure, pass :: set_kbt => set_kbt_xy2d_gpu
     procedure, pass :: set_beta => set_beta_xy2d_gpu
     procedure, pass :: set_initial_magne_autocorrelation_state => set_initial_magne_autocorrelation_state_xy2d_gpu
     procedure, pass, private :: set_autocorrelation_start => set_autocorrelation_start_xy2d_gpu
     !> updater.
     procedure, pass, private :: update_norishiro => update_norishiro_xy2d_gpu
     procedure, pass :: update => update_xy2d_gpu
     procedure, pass :: update_over_relaxation => update_over_relaxation_xy2d_gpu
     !> getter.
     procedure, pass :: nx => nx_xy2d_gpu
     procedure, pass :: ny => ny_xy2d_gpu
     procedure, pass :: nall => nall_xy2d_gpu
     procedure, pass :: kbt => kbt_xy2d_gpu
     procedure, pass :: beta => beta_xy2d_gpu
     procedure, pass :: spins => spins_xy2d_gpu
     !> calculator.
     procedure, pass :: calc_energy_sum => calc_energy_sum_xy2d_gpu
     procedure, pass :: calc_magne_sum => calc_magne_sum_xy2d_gpu
     procedure, pass :: calc_magne_y_sum => calc_magne_y_sum_xy2d_gpu
     ! procedure, pass :: calc_magne_divided_by_initial_magne_sum => calc_magne_divided_by_initial_magne_sum_xy2d_gpu
     procedure, pass :: calc_autocorrelation_sum => calc_autocorrelation_sum_xy2d_gpu
     procedure, pass :: calc_correlation_sum => calc_correlation_sum_xy2d_gpu
  end type xy2d_gpu
contains
  impure subroutine init_xy2d_gpu(this, nx, ny, kbt, iseed)
    class(xy2d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: iseed
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny

    allocate(this%spins_(0 : this%nx_ + 1, 0 : this%ny_ + 1, 1:2))
    allocate(this%autocorrelation_start_spins_(0 : this%nx_ + 1, 0 : this%ny_ + 1, 1:2))
    allocate(this%randoms_(1:this%nx_, 1:this%ny_))
    allocate(this%candidates_(1:this%nx_, 1:this%ny_))
    xy2d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    xy2d_gpu_stat = curandSetPseudoRandomGeneratorSeed(this%rand_gen_, iseed)
    call this%set_allup_spin()
    call this%set_kbt(kbt)
  end subroutine init_xy2d_gpu
 impure subroutine skip_curand_xy2d_gpu(this, n_skip)
    class(xy2d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: n_skip
    if (n_skip == 0_int64) return
    xy2d_gpu_stat = curandSetGeneratorOffset(this%rand_gen_, n_skip)
  end subroutine skip_curand_xy2d_gpu
  !> set_allup_spin_xy2d_gpu: Set 'ferromagnetic' initial state of XY.
  impure subroutine set_allup_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call set_allup_spin_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    call this%update_norishiro()
  end subroutine set_allup_spin_xy2d_gpu
  attributes(global) pure subroutine set_allup_spin_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx
    spins(x, y, :) = [1.0_real64, 0.0_real64] ! 全てX軸方向.
  end subroutine set_allup_spin_sub

  !> set_random_spin_xy2d_gpu: Set 'paramagnetic' initial state of XY.
  !> Set the magnetization vector (M(0), 0) by rotation whole spins.
  impure subroutine set_random_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_, this%randoms_)
    call this%rotate_summation_magne_toward_xaxis()
  end subroutine set_random_spin_xy2d_gpu

  attributes(global) pure subroutine set_random_spin_sub(nall, nx, ny, spins, randoms)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), intent(in) :: randoms(nx, ny)
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx
    spins(x, y, :) = [cos(2 * pi * randoms(x, y)), sin(2 * pi * randoms(x, y))]
  end subroutine set_random_spin_sub

  !> set_finite_magne_spin_xy2d_gpu: Set finite magne initial state of XY.
  !> Set the magnetization vector (M(0), 0) by rotation whole spins.
  impure subroutine set_finite_magne_spin_xy2d_gpu(this, init_magne)
    class(xy2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: init_magne
    real(real64), parameter :: epsilon = 1d-2
    real(real64) :: field_x
    real(real64) :: mx, my, mabs
    call this%set_random_spin()
    field_x = 1d0
    do
       mx = this%calc_magne_sum() / this%nall_
       my = this%calc_magne_y_sum() / this%nall_
       mabs = hypot(mx, my)
       write(error_unit, '(*(g0, 1x))') mabs, init_magne, abs(mabs - init_magne) / init_magne, epsilon
       if (abs(mabs - init_magne) / init_magne < epsilon) exit
       if (mabs > init_magne) then
          field_x = - field_x / 2
       else
          field_x = field_x * 2
       end if
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
       call metropolis_by_field_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
            & (this%nall_, this%nx_, this%ny_, this%spins_, this%randoms_, this%candidates_, field_x, 0d0)
       xy2d_gpu_stat = cudaDeviceSynchronize()
    end do
    call this%rotate_summation_magne_toward_xaxis()
  end subroutine set_finite_magne_spin_xy2d_gpu

  !> set_random_small_spin_xy2d_gpu: Set small 'paramagnetic' initial state of XY.
  !> Set the magnetization vector (M(0), 0) by rotation whole spins.
  impure subroutine set_random_small_spin_xy2d_gpu(this, near_magne)
    class(xy2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: near_magne
    real(real64) :: mx, my, mabs
    call this%set_random_spin()
    do
       mx = this%calc_magne_sum() / this%nall_
       my = this%calc_magne_y_sum() / this%nall_
       mabs = hypot(mx, my)
       write(error_unit, '(*(g0, 1x))') mabs, near_magne
       if (mabs < near_magne) exit
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
       call metropolis_by_field_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
            & (this%nall_, this%nx_, this%ny_, this%spins_, this%randoms_, this%candidates_, -mx, -my)
       xy2d_gpu_stat = cudaDeviceSynchronize()
    end do
    call this%rotate_summation_magne_toward_xaxis()
  end subroutine set_random_small_spin_xy2d_gpu

  !> set_random_near_spin_xy2d_gpu: Set small 'paramagnetic' initial state of XY near `threshold` .
  !> Set the magnetization vector (M(0), 0) by rotation whole spins.
  impure subroutine set_random_near_spin_xy2d_gpu(this, near_magne, diff_parcent)
    class(xy2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: near_magne, diff_parcent
    real(real64) :: mx, my, mabs
    call this%set_random_spin()
    do
       mx = this%calc_magne_sum() / this%nall_
       my = this%calc_magne_y_sum() / this%nall_
       mabs = hypot(mx, my)
       write(error_unit, '(*(g0, 1x))') mabs, near_magne, abs(mabs - near_magne) / near_magne, diff_parcent
       if (abs(mabs - near_magne) / near_magne <= diff_parcent) exit
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
       xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
       call metropolis_by_field_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
            & (this%nall_, this%nx_, this%ny_, this%spins_, this%randoms_, this%candidates_, -mx, -my)
       xy2d_gpu_stat = cudaDeviceSynchronize()
    end do
    call this%rotate_summation_magne_toward_xaxis()
  end subroutine set_random_near_spin_xy2d_gpu

  attributes(global) pure subroutine metropolis_by_field_sub(nall, nx, ny, spins, randoms, candidates, hx, hy)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), intent(in) :: randoms(nx, ny), candidates(nx, ny)
    real(real64), value :: hx, hy
    real(real64) :: candidate(1:2), delta_energy
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx

    candidate(1:2) = [cos(2 * pi * candidates(x, y)), sin(2 * pi * candidates(x, y))]
    delta_energy = - (hx * (candidate(1) - spins(x, y, 1)) + &
         &            hy * (candidate(2) - spins(x, y, 2)))
    if (randoms(x, y) > 1 - exp(delta_energy)) return
    !> randoms(x, y) <= exp(delta_energy)
    spins(x, y, :) = candidate(:)
  end subroutine metropolis_by_field_sub

  !> rotate_summation_magne_toward_xaxis_xy2d_gpu: Rotate whol spins that satisfied \sum_{i} S_{i}^{(y)} = 0.
  impure subroutine rotate_summation_magne_toward_xaxis_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    real(real64) :: mx, my
    real(real64) :: theta
    mx = this%calc_magne_sum()
    my = this%calc_magne_y_sum()
    theta = atan2(my, mx)
    !> (mx, my) == (cos(theta), sin(theta))
    !> Rotate whole spin -theta to set my == 0.
    call rotate_whole_spin_theta_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_, -theta)
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine rotate_summation_magne_toward_xaxis_xy2d_gpu

  !> rotate_summation_magne_and_autocorrelation_toward_xaxis_xy2d_gpu: Rotate whol spins that satisfied \sum_{i} S_{i}^{(y)} = 0.
  impure subroutine rotate_summation_magne_and_autocorrelation_toward_xaxis_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    real(real64) :: mx, my
    real(real64) :: theta
    mx = this%calc_magne_sum()
    my = this%calc_magne_y_sum()
    theta = atan2(my, mx)
    !> (mx, my) == (cos(theta), sin(theta))
    !> Rotate whole spin -theta to set my == 0.
    call rotate_whole_spin_theta_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_, -theta)
    call rotate_whole_spin_theta_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%autocorrelation_start_spins_, -theta)
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine rotate_summation_magne_and_autocorrelation_toward_xaxis_xy2d_gpu

  !> rotate_summation_magne_and_autocorrelation_toward_xaxis_updown_randomly_xy2d_gpu: Rotate whol spins that satisfied \sum_{i} S_{i}^{(y)} = 0.
  impure subroutine rotate_summation_magne_and_autocorrelation_toward_xaxis_updown_randomly_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    real(real64) :: mx, my
    real(real64) :: theta
    mx = this%calc_magne_sum()
    my = this%calc_magne_y_sum()
    theta = atan2(my, mx)

    block
      real(real64), allocatable, device :: r_d(:)
      real(real64) :: r
      allocate(r_d(1))
      xy2d_gpu_stat = curandGenerate(this%rand_gen_, r_d, 1)
      r = r_d(1)
      if (r < 0.5d0) then
         theta = theta + pi
      end if
    end block
    !> (mx, my) == (cos(theta), sin(theta))
    !> Rotate whole spin -theta to set my == 0.
    call rotate_whole_spin_theta_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_, -theta)
    call rotate_whole_spin_theta_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%autocorrelation_start_spins_, -theta)
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine rotate_summation_magne_and_autocorrelation_toward_xaxis_updown_randomly_xy2d_gpu

  attributes(global) pure subroutine rotate_whole_spin_theta_sub(nall, nx, ny, spins, theta)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), value :: theta
    real(real64) :: current_theta
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx
    current_theta = atan2(spins(x, y, 2), spins(x, y, 1))
    spins(x, y, :) = [cos(current_theta + theta), sin(current_theta + theta)]
  end subroutine rotate_whole_spin_theta_sub

  !> update_norishiro_xy2d_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call update_norishiro_updown_sub <<<(this%nx_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    call update_norishiro_leftright_sub <<<(this%ny_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_xy2d_gpu
  attributes(global) pure subroutine update_norishiro_updown_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: x
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (x > nx) return

    ! norishiro top, (1:nx, ny + 1) <- (1:nx, 1)
    spins(x, ny + 1, :) = spins(x, 1, :)
    ! norishiro bottom, (1:nx, 0) <- (1:nx, ny)
    spins(x, 0, :) = spins(x, ny, :)
  end subroutine update_norishiro_updown_sub
  attributes(global) pure subroutine update_norishiro_leftright_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: y
    y = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (y > ny) return
    ! norishiro left, (0, 1:ny) <- (nx, 1:ny)
    spins(0, y, :) = spins(nx, y, :)
    ! norishiro right, (nx + 1, 1:ny) <- (1, 1:ny)
    spins(nx + 1, y, :) = spins(1, y, :)
  end subroutine update_norishiro_leftright_sub

  !> set_kbt_xy2d_gpu: Set new kbt.
  pure subroutine set_kbt_xy2d_gpu(this, kbt)
    class(xy2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_xy2d_gpu
  !> set_beta_xy2d_gpu: Set new kbt.
  pure subroutine set_beta_xy2d_gpu(this, beta)
    class(xy2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
  end subroutine set_beta_xy2d_gpu

  !> set_initial_magne_autocorrelation_state_xy2d_gpu: Set current spin coordination into `autocorrelation_start_spins_`.
  impure subroutine set_initial_magne_autocorrelation_state_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call this%set_autocorrelation_start()
  end subroutine set_initial_magne_autocorrelation_state_xy2d_gpu
  !> set_autocorrelation_start_xy2d_gpu: Set current spin coordination into `autocorrelation_start_spins_`.
  impure subroutine set_autocorrelation_start_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    this%autocorrelation_start_spins_ = this%spins_
  end subroutine set_autocorrelation_start_xy2d_gpu

  !> update_xy2d_gpu: Update by Metropolis method.
  impure subroutine update_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
    call update_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_, this%beta(), this%randoms_, this%candidates_, 0)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call update_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_, this%beta(), this%randoms_, this%candidates_, 1)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_xy2d_gpu
  attributes(global) pure subroutine update_sub(nall, nx, ny, spins, beta, randoms, candidates, offset)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), value :: beta
    real(real64), intent(in) :: randoms(nx, ny), candidates(nx, ny)
    integer(int32), value :: offset
    real(real64) :: candidate(1:2)
    real(real64) :: delta_energy
    integer(int64) :: idx, x, y
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 1
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx + merge(0, 1, iand(y + offset, b'1') == 1)

    candidate(1:2) = [cos(2 * pi * candidates(x, y)), sin(2 * pi * candidates(x, y))]
    delta_energy = calc_delta_energy(nx, ny, spins, x, y, candidate)
    if (randoms(x, y) > exp(- beta * delta_energy)) return
    !> randoms(x, y) <= exp(- beta * delta_energy)
    spins(x, y, :) = candidate(:)
  end subroutine update_sub

  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure real(real64) function calc_delta_energy(nx, ny, spins, x, y, candidate) result(res)
    integer(int64), value :: nx, ny, x, y
    real(real64), intent(in) :: spins(0 : nx + 1, 0 : ny + 1, 1:2), candidate(1:2)
    real(real64) :: center_diff(1:2), neighbor_summ(1:2)
    center_diff = candidate(:) - spins(x, y, :)
    neighbor_summ = spins(x + 1, y, :) + spins(x - 1, y, :) + spins(x, y + 1, :) + spins(x, y - 1, :)
    res = - (center_diff(1) * neighbor_summ(1) + center_diff(2) * neighbor_summ(2))
  end function calc_delta_energy

  !> update_over_relaxation_xy2d_gpu: Update by over relaxation algorithm.
  impure subroutine update_over_relaxation_xy2d_gpu(this, n_steps)
    class(xy2d_gpu), intent(inout) :: this
    integer(int32), intent(in) :: n_steps
    integer(int32) :: i
    do i = 1, n_steps
       call over_relaxation_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_, 0)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call over_relaxation_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_, 1)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
       xy2d_gpu_stat = cudaDeviceSynchronize()
    end do
  end subroutine update_over_relaxation_xy2d_gpu
  !> over_relaxation_sub: Over relaxation method by GPU.
  attributes(global) pure subroutine over_relaxation_sub(nall, nx, ny, spins, offset)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int32), value :: offset
    real(real64) :: local_field(1:2)
    real(real64) :: abs_local_field_inv
    integer(int64) :: idx, x, y
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 1
    if (idx > nall) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx + merge(0, 1, iand(y + offset, b'1') == 1)

    local_field(1:2) = spins(x - 1, y, :) + spins(x + 1, y, :) + spins(x, y - 1, :) + spins(x, y + 1, :)
    abs_local_field_inv = 1 / hypot(local_field(1), local_field(2))
    local_field(1:2) = local_field(1:2) * abs_local_field_inv
    spins(x, y, 1:2) = (2 * sum(local_field(1:2) * spins(x, y, 1:2))) * local_field(1:2) - spins(x, y, 1:2)
    block
      real(real64) :: rabs
      rabs = hypot(spins(x, y, 1), spins(x, y, 2))
      spins(x, y, 1:2) = spins(x, y, 1:2) / rabs
    end block
  end subroutine over_relaxation_sub

  pure integer(int64) function nx_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = this%nx_
  end function nx_xy2d_gpu
  pure integer(int64) function ny_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = this%ny_
  end function ny_xy2d_gpu
  pure integer(int64) function nall_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = this%nall_
  end function nall_xy2d_gpu
  pure real(real64) function kbt_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_xy2d_gpu
  pure real(real64) function beta_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = this%beta_
  end function beta_xy2d_gpu
  pure function spins_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    real(real64), allocatable :: res(:, :, :)
    allocate(res, source = this%spins_)
  end function spins_xy2d_gpu

!> Calculation functions.
  !> calc_energy_sum_xy2d_gpu: Calculate summation of energy.
  impure real(real64) function calc_energy_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_energy_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_)
  end function calc_energy_sum_xy2d_gpu
  !> calc_magne_sum_xy2d_gpu: Calculate summation of magnetization.
  impure real(real64) function calc_magne_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_magne_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_)
  end function calc_magne_sum_xy2d_gpu
  !> calc_magne_y_sum_xy2d_gpu: Calculate summation of y-components of magnetization.
  impure real(real64) function calc_magne_y_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_magne_y_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_)
  end function calc_magne_y_sum_xy2d_gpu
  !> calc_autocorrelation_sum_xy2d_gpu: Calculate autocorrelation between `autocorrelation_start_time` and current time.
  impure real(real64) function calc_autocorrelation_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_autocorrelation_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_, this%autocorrelation_start_spins_)
  end function calc_autocorrelation_sum_xy2d_gpu
  !> calc_correlation_sum_xy2d_gpu: Calculate spin-spin correlation.
  impure real(real64) function calc_correlation_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_correlation_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_)
  end function calc_correlation_sum_xy2d_gpu

!> Implementation of calculation functions.
  !> calc_energy_sum_xy2d_gpu_sub: Calculate the energy density.
  impure real(real64) function calc_energy_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res - spins(x, y, 1) * (spins(x + 1, y, 1) + spins(x, y + 1, 1))
       res = res - spins(x, y, 2) * (spins(x + 1, y, 2) + spins(x, y + 1, 2))
    end do
  end function calc_energy_sum_xy2d_gpu_sub
  !> calc_magne_sum_xy2d_gpu_sub: Calculate summation of x-gradient for magnetization.
  impure real(real64) function calc_magne_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res + spins(x, y, 1)
    end do
  end function calc_magne_sum_xy2d_gpu_sub
  !> calc_magne_y_sum_xy2d_gpu_sub: Calculate summation of y-gradient for magnetization.
  impure real(real64) function calc_magne_y_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res + spins(x, y, 2)
    end do
  end function calc_magne_y_sum_xy2d_gpu_sub
  !> calc_autocorrelation_sum_xy2d_gpu_sub: sub function for `calc_correlation_sum_xy2d_gpu`.
  impure real(real64) function calc_autocorrelation_sum_xy2d_gpu_sub(nall, nx, ny, spins, autocorrelation_start_spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), intent(in), device :: autocorrelation_start_spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res + spins(x, y, 1) * autocorrelation_start_spins(x, y, 1)
       res = res + spins(x, y, 2) * autocorrelation_start_spins(x, y, 2)
    end do
  end function calc_autocorrelation_sum_xy2d_gpu_sub
  !> calc_correlation_sum_xy2d_gpu_sub: sub function for `calc_correlation_sum_xy2d_gpu`.
  impure real(real64) function calc_correlation_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y, next_x, next_y
    res = 0d0
    !$acc parallel loop private(x, y, nx, ny) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       next_y = y + (ny / 2 - 1)
       if (next_y > ny) next_y = next_y - ny
       next_x = x + (nx / 2 - 1)
       if (next_x > nx) next_x = next_x - nx
       res = res + spins(x, y, 1) * spins(next_x, next_y, 1)
       res = res + spins(x, y, 2) * spins(next_x, next_y, 2)
    end do
  end function calc_correlation_sum_xy2d_gpu_sub
end module xy2d_periodic_gpu_m
