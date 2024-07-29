module clock_simple_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  character(len=*), parameter :: version = "GPU_simple"

  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), parameter, public :: mstate = 6
  real(real64), parameter :: pi_state_inv = 2 * pi / mstate

  integer(int64), parameter, public :: nx = 1000_int64, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  real(real64), parameter, public :: kbt = 0.91d0, beta = 1 / kbt

  type(dim3), parameter :: NUM_THREADS = dim3(16_int64, 16_int64, 1_int64)
  type(dim3), parameter :: NUM_BLOCK = dim3(&
       & (nx + NUM_THREADS%x - 1) / NUM_THREADS%x, &
       & (ny + NUM_THREADS%y - 1) / NUM_THREADS%y, &
       & 1_int64)

  integer(int32), allocatable, device :: sixclock(:, :)

  real(real64), allocatable, device :: rnds(:, :, :)
  type(curandGenerator) :: rand_gen

  integer(int64), parameter :: nd = 4

  integer(int32), public, protected :: clock_gpu_stat

  public :: mstate, nx, ny, nall, kbt, beta
  public :: init_sixclock, skip_curand_clock, init_sixclock_order, update_metropolis, calc_energy, calc_magne
contains
  impure subroutine print_version()
    write(output_unit, '(a)') "#"//version
    write(error_unit, '(a)') "#"//version
  end subroutine print_version
 impure subroutine skip_curand_clock(n_skip)
    integer(int64), intent(in) :: n_skip
    if (n_skip == 0_int64) return
    clock_gpu_stat = curandSetGeneratorOffset(rand_gen, n_skip)
  end subroutine skip_curand_clock
  !> init_sixclock: Initialize 6-state clock model, array for random numbers, and random number generator.
  impure subroutine init_sixclock(iseed)
    integer(int32), intent(in) :: iseed
    allocate(sixclock(nx, ny))
    allocate(rnds(1:2, nx, ny))
    clock_gpu_stat = curandCreateGenerator(rand_gen, CURAND_RNG_PSEUDO_XORWOW)
    clock_gpu_stat = curandSetPseudoRandomGeneratorSeed(rand_gen, iseed)
  end subroutine init_sixclock
  !> init_sixclock_order: Set spins with the all-alinged state.
  impure subroutine init_sixclock_order()
    sixclock(:, :) = 0_int32
  end subroutine init_sixclock_order
  !> update_metropolis: Update the lattice with 1MCS.
  impure subroutine update_metropolis()
    clock_gpu_stat = curandGenerate(rand_gen, rnds(:, :, :), 2 * nall)
    !> update even sites.
    call update_sub <<<NUM_BLOCK, NUM_THREADS>>>(sixclock, rnds, 0)
    clock_gpu_stat = cudaDeviceSynchronize()
    !> update odd sites.
    call update_sub <<<NUM_BLOCK, NUM_THREADS>>>(sixclock, rnds, 1)
    clock_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_metropolis
  !> update_sub: Update even or odd sites.
  !> @param sixclock An array in device for spins.
  !> @param rnds An array in device for random numbers.
  !> @param parity_lattice A parity of sites to update. `0` or `1`.
  attributes(global) pure subroutine update_sub(sixclock, rnds, parity_lattice)
    integer(int32), intent(inout) :: sixclock(:, :)
    real(real64), intent(in) :: rnds(:, :, :)
    integer(int32), value :: parity_lattice
    integer(int32) :: new_state
    integer(int64) :: x, y
    integer(int64) :: near_x, near_y
    integer(int32) :: nearest_spins(nd)
    real(real64) :: delta_e
    integer(int32) :: i
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny) return
    if (iand(x + y, b'1') == parity_lattice) return
    !> Calculation about nearest neighbors.
    !> right
    near_x = x + 1
    if (near_x > nx) near_x = 1
    nearest_spins(1) = sixclock(near_x, y)
    !> right
    near_x = x - 1
    if (near_x < 1) near_x = nx
    nearest_spins(2) = sixclock(near_x, y)
    !> up
    near_y = y + 1
    if (near_y > ny) near_y = 1
    nearest_spins(3) = sixclock(x, near_y)
    !> down
    near_y = y - 1
    if (near_y < 1) near_y = ny
    nearest_spins(4) = sixclock(x, near_y)

    !> Calculation of delta_e.
    !> 0 < rnds(1, x, y) <= 1.
    !> 0 < rnds(1, x, y) * (mstate - 1) <= (mstate - -1)
    new_state = sixclock(x, y) + ceiling(rnds(1, x, y) * (mstate - 1))
    if (new_state >= mstate) new_state = new_state - mstate
    delta_e = 0d0
    do i = 1, nd
       delta_e = delta_e + (&
            & - cos((nearest_spins(i) - new_state) * pi_state_inv) &
            & + cos((nearest_spins(i) - sixclock(x, y)) * pi_state_inv))
    end do
    if (delta_e > 0) then
       if (.not. (rnds(2, x, y) <= exp(-beta * delta_e))) return
    end if
    !> delta_e <= 0 .or. rnds(2, x, y) <= exp(-beta * delta_e).
    sixclock(x, y) = new_state
  end subroutine update_sub

  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: i, y, x
    res = 0d0
    !$acc parallel loop private(y, x) reduction(+:res)
    do i = 1, nall
       y = (i - 1) / nx + 1
       x = (i - 1) - (y - 1) * nx + 1 ! == mod(i - 1, nx) + 1
       res = res + cos(sixclock(x, y) * pi_state_inv)
    end do
    res = res * nall_inv
  end function calc_magne
  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: i, y, x, rx, uy
    res = 0d0
    !$acc parallel loop private(y, x, rx, uy) reduction(+:res)
    do i = 1, nall
       y = (i - 1) / nx + 1
       x = (i - 1) - (y - 1) * nx + 1 ! == mod(i - 1, nx) + 1
       rx = x + 1
       if (rx > nx) rx = 1
       uy = y + 1
       if (uy > ny) uy = 1
       res = res - (&
            & cos((sixclock(x, y) - sixclock(rx, y)) * pi_state_inv) + &
            & cos((sixclock(x, y) - sixclock(x, uy)) * pi_state_inv))
    end do
    res = res * nall_inv
  end function calc_energy
end module clock_simple_gpu_m
