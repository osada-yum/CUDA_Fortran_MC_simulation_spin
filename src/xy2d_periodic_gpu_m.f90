!> CUDAFortran implementation for the 2-dimensional XY model with periodic boundary condition.
module xy2d_periodic_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), public, protected :: xy2d_gpu_stat

  type(dim3), parameter :: NUM_THREADS = dim3(32_int64, 16_int64, 1_int64)

  public :: xy2d_gpu
  type :: xy2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     real(real64), allocatable, device :: spins_(:, :, :)
     real(real64), allocatable, device :: randoms_(:, :)
     real(real64), allocatable, device :: candidates_(:, :)
     type(dim3) :: NUM_BLOCK_
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

    this%NUM_BLOCK_ = dim3(&
         & (nx + NUM_THREADS%x - 1) / NUM_THREADS%x, &
         & (ny + NUM_THREADS%y - 1) / NUM_THREADS%y, &
         & 1_int64)

    allocate(this%spins_(0 : this%nx_ + 1, 0:this%ny_ + 1, 1:2))
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
    call set_allup_spin_sub <<<this%NUM_BLOCK_, NUM_THREADS>>>(this%nx_, this%ny_, this%spins_)
    call this%update_norishiro()
  end subroutine set_allup_spin_xy2d_gpu
  attributes(global) pure subroutine set_allup_spin_sub(nx, ny, spins)
    integer(int64), value :: nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: x, y
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return
    spins(x, y, :) = [1.0_real64, 0.0_real64] ! 全てX軸方向.
  end subroutine set_allup_spin_sub

  !> set_random_spin_xy2d_gpu: Set 'paramagnetic' initial state of XY.
  impure subroutine set_random_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub <<<this%NUM_BLOCK_, NUM_THREADS>>>(this%nx_, this%ny_, this%spins_, this%randoms_)
    call this%update_norishiro()
  end subroutine set_random_spin_xy2d_gpu
  attributes(global) pure subroutine set_random_spin_sub(nx, ny, spins, randoms)
    integer(int64), value :: nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), intent(in) :: randoms(nx, ny)
    integer(int64) :: x, y
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return
    spins(x, y, :) = [cos(2 * pi * randoms(x, y)), sin(2 * pi * randoms(x, y))]
  end subroutine set_random_spin_sub

  !> update_norishiro_xy2d_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call update_norishiro_sub <<<this%NUM_BLOCK_, NUM_THREADS>>> &
         & (this%nx_, this%ny_, this%spins_)
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_xy2d_gpu
  attributes(global) pure subroutine update_norishiro_sub(nx, ny, spins)
    integer(int64), value :: nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: x, y
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return

    if (y == 1) then
       ! norishiro top, (1:nx, ny + 1) <- (1:nx, 1)
       spins(x, ny + 1, :) = spins(x, 1, :)
       ! norishiro bottom, (1:nx, 0) <- (1:nx, ny)
       spins(x, 0, :) = spins(x, ny, :)
    else if (x == 1) then
       ! norishiro left, (0, 1:ny) <- (nx, 1:ny)
       spins(0, y, :) = spins(nx, y, :)
       ! norishiro right, (nx + 1, 1:ny) <- (1, 1:ny)
       spins(nx + 1, y, :) = spins(1, y, :)
    end if
  end subroutine update_norishiro_sub

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
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
    call update_sub <<<this%NUM_BLOCK_, NUM_THREADS>>> &
         & (this%nx_, this%ny_, this%spins_, this%beta(), this%randoms_, this%candidates_, 0)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call update_sub <<<this%NUM_BLOCK_, NUM_THREADS>>> &
         & (this%nx_, this%ny_, this%spins_, this%beta(), this%randoms_, this%candidates_, 1)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_xy2d_gpu
  attributes(global) pure subroutine update_sub(nx, ny, spins, beta, randoms, candidates, odd_or_even)
    integer(int64), value :: nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    real(real64), value :: beta
    real(real64), intent(in) :: randoms(nx, ny), candidates(nx, ny)
    integer(int32), value :: odd_or_even
    real(real64) :: candidate(1:2)
    real(real64) :: delta_energy
    integer(int64) :: x, y
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return
    if (iand(x + y, b'1') /= odd_or_even) return

    candidate(1:2) = [cos(2 * pi * candidates(x, y)), sin(2 * pi * candidates(x, y))]
    delta_energy = calc_delta_energy(nx, ny, spins, x, y, candidate)
    if (randoms(x, y) > exp(- beta * delta_energy)) return
    !> randoms(x, y) <= exp(- beta * delta_energy)
    spins(x, y, :) = candidate(:)
  end subroutine update_sub

  !> update_over_relaxation_xy2d_gpu: Update by over relaxation algorithm.
  impure subroutine update_over_relaxation_xy2d_gpu(this, n_steps)
    class(xy2d_gpu), intent(inout) :: this
    integer(int32), intent(in) :: n_steps
    integer(int32) :: i
    do i = 1, n_steps
       call over_relaxation_sub <<<this%NUM_BLOCK_, NUM_THREADS>>> &
         & (this%nx_, this%ny_, this%spins_, 0)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call over_relaxation_sub <<<this%NUM_BLOCK_, NUM_THREADS>>> &
         & (this%nx_, this%ny_, this%spins_, 1)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
       xy2d_gpu_stat = cudaDeviceSynchronize()
    end do
  end subroutine update_over_relaxation_xy2d_gpu
  !> over_relaxation_sub: Over relaxation method by GPU.
  attributes(global) pure subroutine over_relaxation_sub(nx, ny, spins, odd_or_even)
    integer(int64), value :: nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int32), value :: odd_or_even
    real(real64) :: local_field(1:2)
    real(real64) :: abs_local_field_inv
    integer(int64) :: x, y
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return
    if (iand(x + y, b'1') /= odd_or_even) return

    local_field(1:2) = spins(x - 1, y, :) + spins(x + 1, y, :) + spins(x, y - 1, :) + spins(x, y + 1, :)
    abs_local_field_inv = 1 / hypot(local_field(1), local_field(2))
    local_field(1:2) = local_field(1:2) * abs_local_field_inv
    spins(x, y, 1:2) = (2 * sum(local_field(1:2) * spins(x, y, 1:2))) * local_field(1:2) - spins(x, y, 1:2)
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

  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure real(real64) function calc_delta_energy(nx, ny, spins, x, y, candidate) result(res)
    integer(int64), value :: nx, ny, x, y
    real(real64), intent(in) :: spins(0 : nx + 1, 0 : ny + 1, 1:2), candidate(1:2)
    real(real64) :: center_diff(1:2), neighbor_summ(1:2)
    center_diff = candidate(:) - spins(x, y, :)
    neighbor_summ = spins(x + 1, y, :) + spins(x - 1, y, :) + spins(x, y + 1, :) + spins(x, y - 1, :)
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
  pure real(real64) function calc_energy_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_energy_sum_xy2d_gpu_sub(this%nx_, this%ny_, this%spins_)
  end function calc_energy_sum_xy2d_gpu
  pure real(real64) function calc_energy_sum_xy2d_gpu_sub(nx, ny, spins) result(res)
    integer(int64), value :: nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: x, y
    integer(int32) :: i
    res = 0d0
    !$acc parallel loop private(x, y, i) present(spins(:, :, 1:2)) reduction(+:res)
    do i = 1, 2
       do y = 1, ny
          do x = 1, nx
             res = res - spins(x, y, i) * (spins(x + 1, y, i) + spins(x, y + 1, i))
          end do
       end do
    end do
  end function calc_energy_sum_xy2d_gpu_sub
  !> calc_magne_sum_xy2d_gpu: Calculate summation of energy.
  pure real(real64) function calc_magne_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_magne_sum_xy2d_gpu_sub(this%nx_, this%ny_, this%spins_)
  end function calc_magne_sum_xy2d_gpu
  pure real(real64) function calc_magne_sum_xy2d_gpu_sub(nx, ny, spins) result(res)
    integer(int64), value :: nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do y = 1, ny
       do x = 1, nx
          res = res + spins(x, y, 1)
       end do
    end do
  end function calc_magne_sum_xy2d_gpu_sub
end module xy2d_periodic_gpu_m
