!> CUDAFortran implementation for 2-dimensional XY model.
module xy2d_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  use kahan_summation_m
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), public, protected :: xy2d_gpu_stat
  integer(int64), parameter :: NUM_THREADS = 512
  integer(int32), parameter :: lb_exparr = -8, ub_exparr = 8
  public :: xy2d_gpu
  type :: xy2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     integer(int64) :: norishiro_begin_, norishiro_end_
     real(real64), allocatable, device :: spins_(:)
     real(real64), allocatable, device :: randoms_(:)
     real(real64), allocatable, device :: candidates_(:)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_xy2d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_xy2d_gpu
     procedure, pass :: set_random_spin => set_random_spin_xy2d_gpu
     procedure, pass :: set_kbt => set_kbt_xy2d_gpu
     procedure, pass :: set_beta => set_beta_xy2d_gpu
     !> updater.
     procedure, pass, private :: update_norishiro => update_norishiro_xy2d_gpu
     procedure, pass :: update => update_xy2d_gpu
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
  impure subroutine init_xy2d_gpu(this, nx, ny, kbt)
    class(xy2d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%norishiro_begin_ = 1 - nx
    this%norishiro_end_ = this%nall_ + nx
    allocate(this%spins_(this%norishiro_begin_:this%norishiro_end_))
    allocate(this%randoms_(1:this%nall_))
    allocate(this%candidates_(1:this%nall_))
    xy2d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    call this%set_allup_spin()
    call this%set_kbt(kbt)
  end subroutine init_xy2d_gpu
  !> set_allup_spin_xy2d_gpu: Set 'ferromagnetic' initial state of XY.
  pure subroutine set_allup_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    this%spins_(:) = 0.0_real64 !> オールx軸方向.
  end subroutine set_allup_spin_xy2d_gpu
  !> set_random_spin_xy2d_gpu: Set 'paramagnetic' initial state of XY.
  impure subroutine set_random_spin_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub<<<(this%nall_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>>(this%nall_, this%spins_(1:this%nall_), this%randoms_(1:this%nall_))
    call this%update_norishiro()
  end subroutine set_random_spin_xy2d_gpu
  attributes(global) pure subroutine set_random_spin_sub(n, spins, randoms)
    integer(int64), value :: n
    real(real64), intent(inout) :: spins(n)
    real(real64), intent(in) :: randoms(n)
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > n) return !> over norishiro.
    !> 1 <= idx <= this%nall_.
    spins(idx) = 2 * pi * randoms(idx)
  end subroutine set_random_spin_sub
  !> update_norishiro_xy2d_gpu: Update norishiro by GPU.
  pure subroutine update_norishiro_xy2d_gpu(this)
    class(xy2d_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    call update_norishiro_sub <<<(this%nx_ + this%nx_ - 1)/this%nx_, this%nx_>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_)
  end subroutine update_norishiro_xy2d_gpu
  attributes(global) pure subroutine update_norishiro_sub(lb, ub, nx, nall, spins)
    integer(int64), value :: lb, ub, nx, nall
    real(real64), intent(inout) :: spins(lb:ub)
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nx) return !> over norishiro.
    !> 1 <= idx <= this%nx_
    ! norishiro top, [nall+1:nall+nx] <- [1:nx]
    spins(nall + idx) = spins(idx)
    ! norishiro bottom, [-nx:0] <- [nall-nx+1:nall]
    spins(idx - nx) = spins(nall - nx + idx)
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
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_(1:this%nall_), this%nall_)
    xy2d_gpu_stat = curandGenerate(this%rand_gen_, this%candidates_(1:this%nall_), this%nall_)
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%beta(), this%randoms_(1:this%nall_), this%candidates_(1:this%nall_), 1)
    call this%update_norishiro()
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%beta(), this%randoms_(1:this%nall_), this%candidates_(1:this%nall_), 2)
    call this%update_norishiro()
  end subroutine update_xy2d_gpu
  attributes(global) pure subroutine update_sub(lb, ub, nx, nall, spins, beta, randoms, candidates, offset)
    integer(int64), value :: lb, ub, nx, nall
    real(real64), intent(inout) :: spins(lb:ub)
    real(real64), value :: beta
    real(real64), intent(in) :: randoms(nall), candidates(nall)
    integer(int32), value :: offset
    real(real64) :: candidate
    integer(int32) :: delta_energy
    integer(int64) :: idx
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 2 + offset
    if (idx > nall) return
    !> do idx = offset, this%nall_, 2
    candidate = 2 * pi * candidates(idx)
    delta_energy = calc_delta_energy(lb, ub, nx, spins, idx, candidate)
    if (randoms(idx) >= exp(- beta * delta_energy)) return
    !> randoms(idx) < exp(- beta * delta_energy)
    spins(idx) = candidate
  end subroutine update_sub

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
    real(real64), allocatable :: res(:)
    allocate(res, source = this%spins_)
  end function spins_xy2d_gpu

  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure real(real64) function calc_delta_energy(lb, ub, nx, spins, idx, candidate) result(res)
    integer(int64), value :: lb, ub, nx
    real(real64), intent(in) :: spins(lb:ub), candidate
    integer(int64), value :: idx
    res = calc_local_energy(candidate, spins(idx - 1), spins(idx + 1), spins(idx + nx), spins(idx - nx)) &
         & - calc_local_energy(spins(idx), spins(idx - 1), spins(idx + 1), spins(idx + nx), spins(idx - nx))
  contains
    pure real(real64) function calc_local_energy(center, left, right, up, down) result(res)
      real(real64), intent(in) :: center, left, right, up, down
      res = - (cos(center - left) + &
           &   cos(center - right) + &
           &   cos(center - up) + &
           &   cos(center - down))
    end function calc_local_energy
  end function calc_delta_energy
  !> calc_energy_sum_xy2d_gpu: Calculate summation of energy.
  pure real(real64) function calc_energy_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    type(kahan_summation) :: ksum
    real(real64), allocatable :: spins(:)
    integer(int64) :: i
    ksum = kahan_summation(0)
    allocate(spins, source = this%spins_)
    do i = 1, this%nall_
       ksum = ksum + (- cos(spins(i) - spins(i + 1)) + cos(spins(i) - spins(i + this%nx_)))
    end do
    res = ksum%val()
  end function calc_energy_sum_xy2d_gpu
  !> calc_magne_sum_xy2d_gpu: Calculate summation of energy.
  pure real(real64) function calc_magne_sum_xy2d_gpu(this) result(res)
    class(xy2d_gpu), intent(in) :: this
    type(kahan_summation) :: ksum
    real(real64), allocatable :: spins(:)
    integer(int64) :: i
    ksum = kahan_summation(0)
    allocate(spins, source = this%spins_)
    do i = 1, this%nall_
       ksum = ksum + cos(spins(i))
    end do
    res = ksum%val()
  end function calc_magne_sum_xy2d_gpu
end module xy2d_GPU_m
