module clock_dual_lattice_yhalf_tableall_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  character(len=*), parameter :: version = "GPU_dual_lattice_yhalf_tableall"

  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), parameter, public :: mstate = 6
  real(real64), parameter :: pi_state_inv = 2 * pi / mstate

  integer(int64), parameter, public :: nx = 2000_int64, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  real(real64), parameter, public :: kbt = 0.91d0, beta = 1 / kbt

  type(dim3), parameter :: NUM_THREADS = dim3(16_int64, 16_int64, 1_int64)
  type(dim3), parameter :: NUM_BLOCK = dim3(&
       & (nx + NUM_THREADS%x - 1) / NUM_THREADS%x, &
       & (ny / 2 + NUM_THREADS%y - 1) / NUM_THREADS%y, &
       & 1_int64)

  integer(int32), allocatable, device :: sixclock_even(:, :), sixclock_odd(:, :)

  integer(int32) :: global_c, global_r, global_u
  real(real64), parameter :: state_to_magne(0:mstate - 1) = [(cos(global_c * pi_state_inv), global_c = 0, mstate - 1)]
  real(real64), parameter :: state_center_right_up_to_energy(0:mstate - 1, 0:mstate - 1, 0:mstate - 1) = &
       & reshape([ &
       & (&
       &   ( &
       &     (- cos((global_u - global_c) * pi_state_inv) - cos((global_r - global_c) * pi_state_inv), global_c = 0, mstate - 1), &
       & global_u = 0, mstate - 1), &
       & global_r = 0, mstate - 1)], shape = [mstate, mstate, mstate])
  real(real64), allocatable, device :: states_to_prob(:, :, :, :, :, :)

  real(real64), allocatable, device :: rnds(:, :, :)
  type(curandGenerator) :: rand_gen

  integer(int64), parameter :: nd = 4

  integer(int32), public, protected :: clock_gpu_stat

  public :: mstate, nx, ny, nall, kbt, beta
  public :: print_version
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
    integer(int32) :: c, new_c, r, u, l, d
    real(real64) :: delta_e
    real(real64) :: probability_h(0:mstate - 1, 0:mstate - 1, 0:mstate - 1, 0:mstate - 1, 0:mstate - 1, 0:mstate - 1)
    allocate(sixclock_even(nx, ny / 2))
    allocate(sixclock_odd(nx, ny / 2))
    allocate(rnds(1:2, nx, ny))
    clock_gpu_stat = curandCreateGenerator(rand_gen, CURAND_RNG_PSEUDO_XORWOW)
    clock_gpu_stat = curandSetPseudoRandomGeneratorSeed(rand_gen, iseed)
    do d = 0, mstate - 1
       do l = 0, mstate - 1
          do u = 0, mstate - 1
             do r = 0, mstate - 1
                do new_c = 0, mstate - 1
                   do c = 0, mstate - 1
                      delta_e = state_center_right_up_to_energy(new_c, r, u) &
                           & - state_center_right_up_to_energy(c, r, u) &
                           & + state_center_right_up_to_energy(new_c, l, d) &
                           & - state_center_right_up_to_energy(c, l, d)
                      if (delta_e <= 0d0) then
                         probability_h(c, new_c, r, u, l, d) = 1d0
                      else
                         probability_h(c, new_c, r, u, l, d) = exp(- beta * delta_e)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    allocate(states_to_prob, source = probability_h)
  end subroutine init_sixclock
  !> init_sixclock_order: Set spins with the all-alinged state.
  impure subroutine init_sixclock_order()
    sixclock_even(:, :) = 0_int32
    sixclock_odd(:, :) = 0_int32
  end subroutine init_sixclock_order
  !> update_metropolis: Update the lattice with 1MCS.
  impure subroutine update_metropolis()
    clock_gpu_stat = curandGenerate(rand_gen, rnds(:, :, :), 2 * nall)
    !> update even sites.
    call update_sub <<< NUM_BLOCK, NUM_THREADS >>>(sixclock_even, sixclock_odd, rnds, 0)
    clock_gpu_stat = cudaDeviceSynchronize()
    !> update odd sites.
    call update_sub <<< NUM_BLOCK, NUM_THREADS >>>(sixclock_odd, sixclock_even, rnds, 1)
    clock_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_metropolis
  !> update_sub: Update even or odd sites.
  !> @param sixclock_update An array in device for `parity_bit` spins.
  !> @param sixclock_nearest An array in device for `non-parity_bit` spins.
  !> @param rnds An array in device for random numbers.
  !> @param parity_bit A parity of sites to update. `0` or `1`.
  attributes(global) pure subroutine update_sub(sixclock_update, sixclock_nearest, rnds, parity_bit)
    integer(int32), intent(inout) :: sixclock_update(nx, ny / 2)
    integer(int32), intent(in) :: sixclock_nearest(nx, ny / 2)
    real(real64), intent(in) :: rnds(2, nx, ny)
    integer(int32), value :: parity_bit
    integer(int32) :: new_state
    integer(int64) :: x, y
    integer(int64) :: rx, lx, uy, dy, actual_y
    integer(int32) :: nearest_spins(nd)
    integer(int32) :: i
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (x > nx) return
    if (y > ny / 2) return
    !> Calculation about nearest neighbors.
    !> right
    rx = x + 1
    if (rx > nx) rx = 1
    !> left
    lx = x - 1
    if (lx < 1) lx = nx
    !> up
    uy = y + iand(x + 1 + parity_bit, b'1')
    if (uy > ny / 2) uy = 1
    !> down
    dy = y - iand(x + parity_bit, b'1')
    if (dy < 1) dy = ny / 2
    nearest_spins(1) = sixclock_nearest(rx, y)
    nearest_spins(2) = sixclock_nearest(lx, y)
    nearest_spins(3) = sixclock_nearest(x, uy)
    nearest_spins(4) = sixclock_nearest(x, dy)

    !> 0 < rnds(1, x, y) <= 1.
    !> 0 < rnds(1, x, y) * (mstate - 1) <= (mstate - -1)
    actual_y = 2 * y - iand(x + parity_bit, b'1')
    new_state = sixclock_update(x, y) + ceiling(rnds(1, x, actual_y) * (mstate - 1))
    if (new_state >= mstate) new_state = new_state - mstate
    block
      real(real64) :: prob
      prob = states_to_prob(sixclock_update(x, y), new_state, &
           & nearest_spins(1), nearest_spins(3), &
           & nearest_spins(2), nearest_spins(4))
      if (rnds(2, x, actual_y) <= prob) &
           & sixclock_update(x, y) = new_state
    end block
  end subroutine update_sub

  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: y, x
    integer(int64) :: i
    res = calc_magne_sub(sixclock_even, sixclock_odd)
  contains
    pure real(real64) function calc_magne_sub(sixclock_even, sixclock_odd) result(res)
      integer(int32), intent(in), device :: sixclock_even(nx, ny / 2), sixclock_odd(nx, ny / 2)
      integer(int64) :: y, x
      res = 0d0
      !$acc parallel loop private(y, x) reduction(+:res)
      do y = 1, ny / 2
         do x = 1, nx
            res = res + state_to_magne(sixclock_even(x, y)) + state_to_magne(sixclock_odd(x, y))
         end do
      end do
      res = res * nall_inv
    end function calc_magne_sub
  end function calc_magne

  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    res = calc_energy_sub(sixclock_even, sixclock_odd)
  contains
    pure real(real64) function calc_energy_sub(sixclock_even, sixclock_odd) result(res)
      integer(int32), intent(in), device :: sixclock_even(nx, ny / 2), sixclock_odd(nx, ny / 2)
      integer(int64) :: y, x, rx, uy
      res = 0d0
      !$acc parallel loop private(y, x, rx, uy) reduction(+:res)
      do y = 1, ny / 2
         do x = 1, nx
            rx = x + 1
            if (rx > nx) rx = 1
            !> even
            uy = y + iand(x + 1, b'1')
            if (uy > ny / 2) uy = 1
            res = res + state_center_right_up_to_energy(sixclock_even(x, y), sixclock_odd(rx, y), sixclock_odd(x, uy))
            !> odd
            uy = y + iand(x, b'1')
            if (uy > ny / 2) uy = 1
            res = res + state_center_right_up_to_energy(sixclock_odd(x, y), sixclock_even(rx, y), sixclock_even(x, uy))
         end do
      end do
      res = res * nall_inv
    end function calc_energy_sub
  end function calc_energy
end module clock_dual_lattice_yhalf_tableall_gpu_m
