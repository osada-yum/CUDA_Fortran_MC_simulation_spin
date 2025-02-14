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
     ! real(real64), allocatable, device :: autocorrelation_start_spins_(:, :, :)
     ! real(real64), allocatable, device :: randoms_(:, :)
     ! real(real64), allocatable, device :: candidates_(:, :)
     ! type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_xy2d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_xy2d_gpu
     ! procedure, pass :: set_initial_magne_autocorrelation_state => set_initial_magne_autocorrelation_state_xy2d_gpu
     ! procedure, pass, private :: set_autocorrelation_start => set_autocorrelation_start_xy2d_gpu
     !> getter.
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

    allocate(this%spins_(0 : this%nx_ + 1, 0 : this%ny_ + 1, 1:2))
    ! allocate(this%autocorrelation_start_spins_(0 : this%nx_ + 1, 0 : this%ny_ + 1, 1:2))
    ! allocate(this%randoms_(1:this%nx_, 1:this%ny_))
    ! allocate(this%candidates_(1:this%nx_, 1:this%ny_))
    ! xy2d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    ! xy2d_gpu_stat = curandSetPseudoRandomGeneratorSeed(this%rand_gen_, iseed)
    call this%set_allup_spin()
  end subroutine init_xy2d_gpu
  !> set_allup_spin_xy2d_gpu: Set 'ferromagnetic' initial state of XY.
  impure subroutine set_allup_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    this%spins_(:, :, :) = 1d0
    ! call set_allup_spin_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
    !      & (this%nall_, this%nx_, this%ny_, this%spins_)
  end subroutine set_allup_spin_xy2d_gpu
  attributes(global) impure subroutine set_allup_spin_sub(nall, nx, ny)
    integer(int64), value :: nall, nx, ny
    integer(int64) :: idx, x, y
  end subroutine set_allup_spin_sub

  attributes(global) impure subroutine update_norishiro_updown_sub(nall, nx, ny)
    integer(int64), value :: nall, nx, ny
    integer(int64) :: x
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
  end subroutine update_norishiro_updown_sub
  attributes(global) impure subroutine update_norishiro_leftright_sub(nall, nx, ny)
    integer(int64), value :: nall, nx, ny
    integer(int64) :: y
    y = (blockIdx%x - 1) * blockDim%x + threadIdx%x
  end subroutine update_norishiro_leftright_sub

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
    !$acc parallel loop private(x, y) present(spins) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       ! x = idx - (y - 1) * nx
       res = res + y
       ! res = res - spins(x, y, 1) * (spins(x + 1, y, 1) + spins(x, y + 1, 1))
       ! res = res - spins(x, y, 2) * (spins(x + 1, y, 2) + spins(x, y + 1, 2))
    end do
  end function calc_energy_sum_xy2d_gpu_sub
  !> calc_magne_sum_xy2d_gpu_sub: Calculate summation of x-gradient for magnetization.
  impure real(real64) function calc_magne_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins) reduction(+:res)
    do idx = 1, nall
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res + spins(x, y, 1)
    end do
  end function calc_magne_sum_xy2d_gpu_sub
  !> calc_correlation_sum_xy2d_gpu_sub: sub function for `calc_correlation_sum_xy2d_gpu`.
  impure real(real64) function calc_correlation_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), intent(in) :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny + 1, 1:2)
    integer(int64) :: idx, x, y, next_x, next_y
    res = 0d0
    !$acc parallel loop private(x, y, nx, ny) present(spins) reduction(+:res)
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
