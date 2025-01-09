!> CUDAFortran implementation for the 2-dimensional XY model with periodic boundary condition.
!> Spins are represented by two arrays.
!> Split spins into even and odd for speedup.
module xy2d_periodic_yhalf_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), public, protected :: xy2d_gpu_stat

  integer(int32), parameter :: nd = 4
  !> up, down, right, left
  integer(int64), parameter :: dy(nd, 0:1) = &
       & reshape([[1_int64, 0_int64, 0_int64, 0_int64], [0_int64, -1_int64, 0_int64, 0_int64]], shape = [4, 2])
  integer(int64), parameter :: dx(nd) = [0_int64, 0_int64, 1_int64, -1_int64]

  integer(int64), parameter :: NUM_THREADS = 512

  public :: xy2d_gpu
  type :: xy2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     real(real64), allocatable, device :: spins_even_(:, :, :), spins_odd_(:, :, :)
     real(real64), allocatable, device :: randoms_(:, :)
     real(real64), allocatable, device :: candidates_(:, :)
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

    allocate(this%spins_even_(0 : this%nx_ + 1, 0:this%ny_ / 2 + 1, 1:2))
    allocate(this%spins_odd_(0 : this%nx_ + 1, 0:this%ny_ / 2 + 1, 1:2))
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
    call set_allup_spin_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_even_)
    call set_allup_spin_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_)
    call this%update_norishiro()
  end subroutine set_allup_spin_xy2d_gpu
  attributes(global) pure subroutine set_allup_spin_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall / 2) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx
    spins(x, y, :) = [1.0_real64, 0.0_real64] ! 全てX軸方向.
  end subroutine set_allup_spin_sub

  !> set_random_spin_xy2d_gpu: Set 'paramagnetic' initial state of XY.
  impure subroutine set_random_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_even_, this%randoms_(:, 1 : this%ny_ / 2))
    call set_random_spin_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_, this%randoms_(:, this%ny_ / 2 + 1 : this%ny_))
    call this%update_norishiro()
  end subroutine set_random_spin_xy2d_gpu
  attributes(global) pure subroutine set_random_spin_sub(nall, nx, ny, spins, randoms)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), intent(in) :: randoms(1:nx, 1:ny/2)
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall / 2) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx
    spins(x, y, :) = [cos(2 * pi * randoms(x, y)), sin(2 * pi * randoms(x, y))]
  end subroutine set_random_spin_sub

  !> update_norishiro_xy2d_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    call update_norishiro_updown_sub <<<(this%nx_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_even_)
    call update_norishiro_updown_sub <<<(this%nx_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_)
    call update_norishiro_leftright_sub <<<(this%ny_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_even_)
    call update_norishiro_leftright_sub <<<(this%ny_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_)
    xy2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_xy2d_gpu
  attributes(global) pure subroutine update_norishiro_updown_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int64) :: x
    x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (x > nx) return
    ! norishiro bottom, (1:nx, 0) <- (1:nx, ny/2)
    spins(x, 0, 1) = spins(x, ny / 2, 1)
    spins(x, 0, 2) = spins(x, ny / 2, 2)
    ! norishiro top, (1:nx, ny/2 + 1) <- (1:nx, 1)
    spins(x, ny / 2 + 1, 1) = spins(x, 1, 1)
    spins(x, ny / 2 + 1, 2) = spins(x, 1, 2)
  end subroutine update_norishiro_updown_sub
  attributes(global) pure subroutine update_norishiro_leftright_sub(nall, nx, ny, spins)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int64) :: y
    y = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (y > ny / 2) return
    ! norishiro left, (0, 1:ny/2) <- (nx, 1:ny/2)
    spins(0, y, 1:2) = spins(nx, y, 1:2)
    ! norishiro right, (nx + 1, 1:ny/2) <- (1, 1:ny/2)
    spins(nx + 1, y, 1:2) = spins(1, y, 1:2)
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
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_, this%nall_)
    call update_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_even_, this%spins_odd_, this%beta(), &
         & this%randoms_(:, 1:this%ny_ / 2), this%candidates_(:, 1 : this%ny_ / 2), 0)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    call update_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_, this%spins_even_, this%beta(), &
         & this%randoms_(:, this%ny_ / 2 + 1 : this%ny_), this%candidates_(:, this%ny_ / 2 + 1 : this%ny_), 1)
    xy2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine update_xy2d_gpu
  attributes(global) pure subroutine update_sub(nall, nx, ny, spins_update, spins_around, beta, randoms, candidates, even_or_odd)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins_update(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), intent(in) :: spins_around(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), value :: beta
    real(real64), intent(in) :: randoms(1:nx, 1:ny/2), candidates(1:nx, 1:ny/2)
    integer(int32), value :: even_or_odd !> 0 if even, 1 if odd.
    real(real64) :: candidate(1:2)
    real(real64) :: delta_energy
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall / 2) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx

    call sincospi(2 * candidates(x, y), candidate(2), candidate(1))
    delta_energy = calc_delta_energy(nx, ny, spins_update, spins_around, x, y, candidate, even_or_odd)
    if (randoms(x, y) > exp(- beta * delta_energy)) return
    !> randoms(x, y) <= exp(- beta * delta_energy)
    spins_update(x, y, :) = candidate(:)
  end subroutine update_sub

  !> update_over_relaxation_xy2d_gpu: Update by over relaxation algorithm.
  impure subroutine update_over_relaxation_xy2d_gpu(this, n_steps)
    class(xy2d_gpu), intent(inout) :: this
    integer(int32), intent(in) :: n_steps
    integer(int32) :: i
    do i = 1, n_steps
       call over_relaxation_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_even_, this%spins_odd_, 0)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
       call over_relaxation_sub <<<(this%nall_ / 2 + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (this%nall_, this%nx_, this%ny_, this%spins_odd_, this%spins_even_, 1)
       xy2d_gpu_stat = cudaDeviceSynchronize()
       call this%update_norishiro()
    end do
  end subroutine update_over_relaxation_xy2d_gpu
  !> over_relaxation_sub: Over relaxation method by GPU.
  attributes(global) pure subroutine over_relaxation_sub(nall, nx, ny, spins_update, spins_around, even_or_odd)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(inout) :: spins_update(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), intent(in) :: spins_around(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int32), value :: even_or_odd
    real(real64) :: local_field(1:2)
    real(real64) :: abs_local_field_inv
    integer(int64) :: idx, x, y
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nall / 2) return
    y = (idx - 1) / nx + 1
    x = idx - (y - 1) * nx

    local_field(1:2) = calc_local_field(nx, ny, spins_around, x, y, even_or_odd)
    abs_local_field_inv = 1 / hypot(local_field(1), local_field(2))
    local_field(1:2) = local_field(1:2) * abs_local_field_inv
    spins_update(x, y, 1:2) = (2 * sum(local_field(1:2) * spins_update(x, y, 1:2))) * local_field(1:2) - spins_update(x, y, 1:2)
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

  !> calc_local_field  Cauculate local field that is the summation of spin variables around the center spin.
  attributes(device) pure function calc_local_field(nx, ny, spins_around, x, y, even_or_odd) result(res)
    real(real64) :: res(1:2)
    integer(int64), value :: nx, ny, x, y
    real(real64), intent(in) :: spins_around(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int32), intent(in) :: even_or_odd
    integer(int64) :: near_x, near_y
    integer(int32) :: d
    res(1:2) = 0.0d0
    do d = 1, nd / 2
       near_y = y + dy(d, iand(x + even_or_odd, b'1'))
       res(1:2) = res(1:2) + spins_around(x, near_y, 1:2)
    end do
    do d = nd / 2 + 1, nd
       near_x = x + dx(d)
       res(1:2) = res(1:2) + spins_around(near_x, y, 1:2)
    end do
  end function calc_local_field
  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure real(real64) function calc_delta_energy(nx, ny, spins_update, spins_around, x, y, candidate, even_or_odd) result(res)
    integer(int64), value :: nx, ny, x, y
    real(real64), intent(in) :: spins_update(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), intent(in) :: spins_around(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    real(real64), intent(in) :: candidate(1:2)
    integer(int32), intent(in) :: even_or_odd
    real(real64) :: center_diff(1:2), neighbor_summ(1:2)
    center_diff(1:2) = candidate(:) - spins_update(x, y, :)
    neighbor_summ(1:2) = calc_local_field(nx, ny, spins_around, x, y, even_or_odd)
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
    res = calc_energy_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_even_) + &
         & calc_energy_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_odd_)
  end function calc_energy_sum_xy2d_gpu
  pure real(real64) function calc_energy_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall / 2
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res - spins(x, y, 1) * (spins(x + 1, y, 1) + spins(x, y + 1, 1))
       res = res - spins(x, y, 2) * (spins(x + 1, y, 2) + spins(x, y + 1, 2))
    end do
  end function calc_energy_sum_xy2d_gpu_sub
  !> calc_magne_sum_xy2d_gpu: Calculate summation of energy.
  pure real(real64) function calc_magne_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    res = calc_magne_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_even_) + &
         & calc_magne_sum_xy2d_gpu_sub(this%nall_, this%nx_, this%ny_, this%spins_odd_)
  end function calc_magne_sum_xy2d_gpu
  pure real(real64) function calc_magne_sum_xy2d_gpu_sub(nall, nx, ny, spins) result(res)
    integer(int64), value :: nall, nx, ny
    real(real64), intent(in), device :: spins(0 : nx + 1, 0 : ny / 2 + 1, 1:2)
    integer(int64) :: idx, x, y
    res = 0d0
    !$acc parallel loop private(x, y) present(spins(:, :, 1:2)) reduction(+:res)
    do idx = 1, nall / 2
       y = (idx - 1) / nx + 1
       x = idx - (y - 1) * nx
       res = res + spins(x, y, 1)
    end do
  end function calc_magne_sum_xy2d_gpu_sub
end module xy2d_periodic_yhalf_gpu_m
