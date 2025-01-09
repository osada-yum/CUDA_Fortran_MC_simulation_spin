!> Wrong calculation?
module xy2d_periodic_samples_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), public, protected :: xy2d_gpu_stat

  integer(int64), parameter, public :: NUM_THREADS = 32

  public :: xy2d_gpu
  type :: xy2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     real(real64), allocatable, device :: spins_(:, :, :, :)
     real(real64), allocatable, device :: randoms_(:, :, :)
     real(real64), allocatable, device :: candidates_(:, :, :)
     real(real64), allocatable, device :: summ_ene_(:), summ_mag_(:)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_xy2d_gpu
     procedure, pass :: skip_curand => skip_curand_xy2d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_xy2d_gpu
     procedure, pass :: set_random_spin => set_random_spin_xy2d_gpu
     procedure, pass :: set_kbt => set_kbt_xy2d_gpu
     procedure, pass :: set_beta => set_beta_xy2d_gpu
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

    allocate(this%spins_(1:2, NUM_THREADS, 0 : this%nx_ + 1, 0:this%ny_ + 1))
    allocate(this%randoms_(NUM_THREADS, 1:this%nx_, 1:this%ny_))
    allocate(this%candidates_(NUM_THREADS, 1:this%nx_, 1:this%ny_))
    xy2d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    xy2d_gpu_stat = curandSetPseudoRandomGeneratorSeed(this%rand_gen_, iseed)

    allocate(this%summ_ene_(NUM_THREADS), source = 0d0)
    allocate(this%summ_mag_(NUM_THREADS), source = 0d0)

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
    call set_allup_spin_sub <<<1, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    call this%update_norishiro()
  end subroutine set_allup_spin_xy2d_gpu
  attributes(global) pure subroutine set_allup_spin_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    integer(int64) :: sample, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do y = 1, ny
       do x = 1, nx
          spins(:, sample, x, y) = [1.0_real64, 0.0_real64] ! 全てX軸方向.
       end do
    end do
  end subroutine set_allup_spin_sub

  !> set_random_spin_xy2d_gpu: Set 'paramagnetic' initial state of XY.
  impure subroutine set_random_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_ * NUM_THREADS)
    call set_random_spin_sub <<<1, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_, this%randoms_)
    call this%update_norishiro()
  end subroutine set_random_spin_xy2d_gpu
  attributes(global) pure subroutine set_random_spin_sub(nall, nx, ny, spins, randoms)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    real(real64), intent(in) :: randoms(1:NUM_THREADS, nx, ny)
    integer(int64) :: sample, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do y = 1, ny
       do x = 1, nx
          spins(:, sample, x, y) = [cos(2 * pi * randoms(sample, x, y)), sin(2 * pi * randoms(sample, x, y))]
       end do
    end do
  end subroutine set_random_spin_sub

  !> update_norishiro_xy2d_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call update_norishiro_updown_sub <<<1, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    call update_norishiro_leftright_sub <<<1, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_)
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_xy2d_gpu
  attributes(global) pure subroutine update_norishiro_updown_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    integer(int64) :: sample, x
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do x = 1, nx
       ! norishiro top, (1:nx, ny + 1) <- (1:nx, 1)
       spins(:, sample, x, ny + 1) = spins(:, sample, x, 1)
       ! norishiro bottom, (1:nx, 0) <- (1:nx, ny)
       spins(:, sample, x, 0) = spins(:, sample, x, ny)
    end do
  end subroutine update_norishiro_updown_sub
  attributes(global) pure subroutine update_norishiro_leftright_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    integer(int64) :: sample, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do y = 1, ny
       ! norishiro left, (0, 1:ny) <- (nx, 1:ny)
       spins(:, sample, 0, y) = spins(:, sample, nx, y)
       ! norishiro right, (nx + 1, 1:ny) <- (1, 1:ny)
       spins(:, sample, nx + 1, y) = spins(:, sample, 1, y)
    end do
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

  !> update_xy2d_gpu: Update by Metropolis method.
  impure subroutine update_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_ * NUM_THREADS)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_ * NUM_THREADS)

    call update_sub <<<1, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_, this%beta(), this%randoms_, this%candidates_)
    xy2d_gpu_stat = cudaDeviceSynchronize()

    call this%update_norishiro()
  end subroutine update_xy2d_gpu
  attributes(global) pure subroutine update_sub(nall, nx, ny, spins, beta, randoms, candidates)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    real(real64), value :: beta
    real(real64), intent(in) :: randoms(1:NUM_THREADS, nx, ny), candidates(1:NUM_THREADS, nx, ny)
    real(real64) :: candidate(1:2)
    real(real64) :: delta_energy
    integer(int64) :: sample, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do y = 1, ny
       do x = 2 - iand(y, b'1'), nx, 2
          candidate(1:2) = [cos(2 * pi * candidates(sample, x, y)), sin(2 * pi * candidates(sample, x, y))]
          delta_energy = calc_delta_energy(nx, ny, sample, spins, x, y, candidate)
          if (randoms(sample, x, y) > exp(- beta * delta_energy)) cycle
          !> randoms(x, y, sample) <= exp(- beta * delta_energy)n
          spins(:, sample, x, y) = candidate(:)
       end do
    end do
  end subroutine update_sub

  !> update_over_relaxation_xy2d_gpu: Update by over relaxation algorithm.
  impure subroutine update_over_relaxation_xy2d_gpu(this, n_steps)
    class(xy2d_gpu), intent(inout) :: this
    integer(int32), intent(in) :: n_steps
    integer(int32) :: i
    do i = 1, n_steps
       call over_relaxation_sub <<<1, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
    end do
  end subroutine update_over_relaxation_xy2d_gpu
  !> over_relaxation_sub: Over relaxation method by GPU.
  attributes(global) pure subroutine over_relaxation_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    real(real64) :: local_field(1:2)
    real(real64) :: abs_local_field_inv
    integer(int64) :: sample, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    do y = 1, ny
       do x = 2 - iand(y, b'1'), nx, 2
          local_field(1:2) = spins(1:2, sample, x - 1, y) + spins(1:2, sample, x + 1, y) + spins(1:2, sample, x, y - 1) + spins(1:2, sample, x, y + 1)
          abs_local_field_inv = 1 / hypot(local_field(1), local_field(2))
          local_field(1:2) = local_field(1:2) * abs_local_field_inv
          spins(1:2, sample, x, y) = (2 * sum(local_field(1:2) * spins(1:2, sample, x, y))) * local_field(1:2) - spins(1:2, sample, x, y)
       end do
    end do
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
    real(real64), allocatable :: res(:, :, :, :)
    allocate(res, source = this%spins_)
  end function spins_xy2d_gpu

  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure real(real64) function calc_delta_energy(nx, ny, sample, spins, x, y, candidate) result(res)
    integer(int64), value :: nx, ny, sample, x, y
    real(real64), intent(in) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1), candidate(1:2)
    real(real64) :: center_diff(1:2), neighbor_summ(1:2)
    center_diff = candidate(:) - spins(1:2, sample, x, y)
    neighbor_summ = spins(1:2, sample, x + 1, y) + spins(1:2, sample, x - 1, y) + spins(1:2, sample, x, y + 1) + spins(1:2, sample, x, y - 1)
    res = - (center_diff(1) * neighbor_summ(1) + center_diff(2) * neighbor_summ(2))
  end function calc_delta_energy
  attributes(device) pure real(real64) function calc_local_energy(center, left, right, up, down) result(res)
    real(real64), value :: center, left, right, up, down
    res = - (cos(center - left) + &
         &   cos(center - right) + &
         &   cos(center - up) + &
         &   cos(center - down))
  end function calc_local_energy

  !> calc_energy_sum_xy2d_gpu: Calculate summation of energy.
  subroutine calc_energy_sum_xy2d_gpu(this, res)
    class(xy2d_gpu), intent(in) :: this
    real(real64), intent(inout) :: res(NUM_THREADS)
    call calc_energy_sum_xy2d_gpu_sub <<<1, NUM_THREADS>>> (this%nall_, this%nx_, this%ny_, this%spins_, this%summ_ene_)
    res = this%summ_ene_
  end subroutine calc_energy_sum_xy2d_gpu
  attributes(global) subroutine calc_energy_sum_xy2d_gpu_sub(nall, nx, ny, spins, summ)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(in) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    real(real64), intent(inout) :: summ(1:NUM_THREADS)
    integer(int64) :: sample, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (sample > NUM_THREADS) return

    summ(sample) = 0d0
    do y = 1, ny
       do x = 1, nx
          summ(sample) = summ(sample) - spins(1, sample, x, y) * (spins(1, sample, x + 1, y) + spins(1, sample, x, y + 1))
          summ(sample) = summ(sample) - spins(2, sample, x, y) * (spins(2, sample, x + 1, y) + spins(2, sample, x, y + 1))
       end do
    end do
  end subroutine calc_energy_sum_xy2d_gpu_sub

  !> calc_magne_sum_xy2d_gpu: Calculate summation of energy.
  subroutine calc_magne_sum_xy2d_gpu(this, res)
    class(xy2d_gpu), intent(in) :: this
    real(real64), intent(inout) :: res(NUM_THREADS)
    call calc_magne_sum_xy2d_gpu_sub <<<1, NUM_THREADS>>> (this%nall_, this%nx_, this%ny_, this%spins_, this%summ_mag_)
    res = this%summ_mag_
  end subroutine calc_magne_sum_xy2d_gpu
  attributes(global) subroutine calc_magne_sum_xy2d_gpu_sub(nall, nx, ny, spins, summ)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(in) :: spins(1:2, 1:NUM_THREADS, 0 : nx + 1, 0 : ny + 1)
    real(real64), intent(inout) :: summ(1:NUM_THREADS)
    integer(int64) :: sample, pos, x, y
    sample = (blockIdx%x - 1) * blockDim%x + threadIdx%x

    summ(sample) = 0d0
    do y = 1, ny
       do x = 1, nx
          summ(sample) = summ(sample) + spins(1, sample, x, y)
       end do
    end do
  end subroutine calc_magne_sum_xy2d_gpu_sub
end module xy2d_periodic_samples_gpu_m
