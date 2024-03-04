!> CUDAFortran implementation for 2-dimensional Ising model.
module ising2d_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  integer(int32), public, protected :: ising2d_gpu_stat
  integer(int64), parameter :: NUM_THREADS = 512
  integer(int32), parameter :: lb_exparr = -8, ub_exparr = 8
  public :: ising2d_gpu
  type :: ising2d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     integer(int64) :: norishiro_begin_, norishiro_end_
     integer(int32), allocatable, device :: spins_(:)
     real(real64), allocatable, device :: exparr_(:)
     real(real64), allocatable, device :: randoms_(:)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_ising2d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_ising2d_gpu
     procedure, pass :: set_random_spin => set_random_spin_ising2d_gpu
     procedure, pass :: set_kbt => set_kbt_ising2d_gpu
     procedure, pass :: set_beta => set_beta_ising2d_gpu
     !> updater.
     procedure, pass, private :: update_exparr => update_exparr_ising2d_gpu
     procedure, pass, private :: update_norishiro => update_norishiro_ising2d_gpu
     procedure, pass :: update => update_ising2d_gpu
     !> getter.
     procedure, pass :: nx => nx_ising2d_gpu
     procedure, pass :: ny => ny_ising2d_gpu
     procedure, pass :: nall => nall_ising2d_gpu
     procedure, pass :: kbt => kbt_ising2d_gpu
     procedure, pass :: beta => beta_ising2d_gpu
     procedure, pass :: spins => spins_ising2d_gpu
     !> calculator.
     procedure, pass :: calc_energy_sum => calc_energy_sum_ising2d_gpu
     procedure, pass :: calc_magne_sum => calc_magne_sum_ising2d_gpu
  end type ising2d_gpu
contains
  impure subroutine init_ising2d_gpu(this, nx, ny, kbt)
    class(ising2d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%norishiro_begin_ = 1 - nx
    this%norishiro_end_ = this%nall_ + nx
    allocate(this%spins_(this%norishiro_begin_:this%norishiro_end_))
    allocate(this%randoms_(1:this%nall_))
    ising2d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    call this%set_allup_spin()
    allocate(this%exparr_(lb_exparr:ub_exparr))
    call this%set_kbt(kbt)
  end subroutine init_ising2d_gpu
  !> set_allup_spin_ising2d_gpu: Set 'ferromagnetic' initial state of Ising.
  pure subroutine set_allup_spin_ising2d_gpu(this)
    class(ising2d_gpu), intent(inout) :: this
    this%spins_(:) = 1_int32
  end subroutine set_allup_spin_ising2d_gpu
  !> set_random_spin_ising2d_gpu: Set 'paramagnetic' initial state of Ising.
  impure subroutine set_random_spin_ising2d_gpu(this)
    class(ising2d_gpu), intent(inout) :: this
    ising2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub<<<(this%nall_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>>(this%nall_, this%spins_(1:this%nall_), this%randoms_(1:this%nall_))
    ising2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine set_random_spin_ising2d_gpu
  attributes(global) pure subroutine set_random_spin_sub(n, spins, randoms)
    integer(int64), value :: n
    integer(int32), intent(inout) :: spins(n)
    real(real64), intent(in) :: randoms(n)
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > n) return !> over norishiro.
    !> 1 <= idx <= this%nall_.
    spins(idx) = merge(1, -1, randoms(idx) < 0.5_real64)
  end subroutine set_random_spin_sub
  !> update_norishiro_ising2d_gpu: Update norishiro by GPU.
  pure subroutine update_norishiro_ising2d_gpu(this)
    class(ising2d_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    call update_norishiro_sub <<<(this%nx_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_)
    ising2d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_ising2d_gpu
  attributes(global) pure subroutine update_norishiro_sub(lb, ub, nx, nall, spins)
    integer(int64), value :: lb, ub, nx, nall
    integer(int32), intent(inout) :: spins(lb:ub)
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > nx) return !> over norishiro.
    !> 1 <= idx <= this%nx_
    ! norishiro top, [nall+1:nall+nx] <- [1:nx]
    spins(nall + idx) = spins(idx)
    ! norishiro bottom, [-nx:0] <- [nall-nx+1:nall]
    spins(idx - nx) = spins(nall - nx + idx)
  end subroutine update_norishiro_sub
  !> set_kbt_ising2d_gpu: Set new kbt.
  pure subroutine set_kbt_ising2d_gpu(this, kbt)
    class(ising2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_ising2d_gpu
  !> set_beta_ising2d_gpu: Set new kbt.
  pure subroutine set_beta_ising2d_gpu(this, beta)
    class(ising2d_gpu), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
    call this%update_exparr()
  end subroutine set_beta_ising2d_gpu
  !> update_exparr_ising2d_gpu: Update exparr_(:) by temperature.
  !> Because Ising model is discrete, we can calculate `exp(-βΔE)` in advance.
  pure subroutine update_exparr_ising2d_gpu(this)
    class(ising2d_gpu), intent(inout) :: this
    real(real64) :: tmp
    integer(int32) :: diff
    this%exparr_(:) = 1.0_real64
    do diff = 1, ub_exparr
       !> diff > 0.
       this%exparr_(diff) = exp(- this%beta() * diff)
    end do
  end subroutine update_exparr_ising2d_gpu
  !> update_ising2d_gpu: Update by Metropolis method.
  impure subroutine update_ising2d_gpu(this)
    class(ising2d_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    ising2d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_(1:this%nall_), this%nall_)
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%randoms_(1:this%nall_), this%exparr_, 1)
    ising2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%randoms_(1:this%nall_), this%exparr_, 2)
    ising2d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine update_ising2d_gpu
  attributes(global) pure subroutine update_sub(lb, ub, nx, nall, spins, randoms, exparr, offset)
    integer(int64), value :: lb, ub, nx, nall
    integer(int32), intent(inout) :: spins(lb:ub)
    real(real64), intent(in) :: randoms(nall), exparr(lb_exparr:ub_exparr)
    integer(int32), value :: offset
    integer(int32) :: delta_energy
    integer(int64) :: idx
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 2 + offset
    if (idx > nall) return
    !> do idx = offset, this%nall_, 2
    delta_energy = calc_delta_energy(lb, ub, nx, spins, idx)
    if (randoms(idx) >= exparr(delta_energy)) return
    !> randoms(idx) < exparr(delta_energy)
    spins(idx) = - spins(idx)
  end subroutine update_sub

  pure integer(int64) function nx_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = this%nx_
  end function nx_ising2d_gpu
  pure integer(int64) function ny_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = this%ny_
  end function ny_ising2d_gpu
  pure integer(int64) function nall_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = this%nall_
  end function nall_ising2d_gpu
  pure real(real64) function kbt_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_ising2d_gpu
  pure real(real64) function beta_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = this%beta_
  end function beta_ising2d_gpu
  pure function spins_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    integer(int32), allocatable :: res(:)
    allocate(res, source = this%spins_)
  end function spins_ising2d_gpu

  !> calc_delta_energy: Calculate delta energy if spins_(idx) is flipped.
  attributes(device) pure integer(int32) function calc_delta_energy(lb, ub, nx, spins, idx) result(res)
    integer(int64), value :: lb, ub, nx
    integer(int32), intent(in) :: spins(lb:ub)
    integer(int64), value :: idx
    res = 2 * spins(idx) * (spins(idx + 1) + spins(idx - 1) + spins(idx + nx) + spins(idx - nx))
  end function calc_delta_energy
  !> calc_energy_sum_ising2d_gpu: Calculate summation of energy.
  pure integer(int64) function calc_energy_sum_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = calc_energy_sum_sub(this%nall_, this%norishiro_end_, this%spins_(1:this%norishiro_end_))
  contains
    pure integer(int64) function calc_energy_sum_sub(n, norishiro_end, spins) result(res)
      integer(int64), intent(in) :: n, norishiro_end
      integer(int32), intent(in), device :: spins(norishiro_end)
      integer(int64) :: i
      res = 0_int64
      !$acc parallel loop present(spins) reduction(+:res)
      do i = 1, n
         res = res - int(this%spins_(i) * (this%spins_(i + 1) + this%spins_(i + this%nx_)), int64)
      end do
    end function calc_energy_sum_sub
  end function calc_energy_sum_ising2d_gpu
  !> calc_magne_sum_ising2d_gpu: Calculate summation of energy.
  pure integer(int64) function calc_magne_sum_ising2d_gpu(this) result(res)
    class(ising2d_gpu), intent(in) :: this
    res = calc_magne_sum_sub(this%nall_, this%spins_(1:this%nall_))
  contains
    pure integer(int64) function calc_magne_sum_sub(n, spins) result(res)
      integer(int64), intent(in) :: n
      integer(int32), intent(in), device :: spins(n)
      integer(int64) :: i
      res = 0_int64
      !$acc parallel loop present(spins) reduction(+:res)
      do i = 1, n
         res = res + spins(i)
      end do
    end function calc_magne_sum_sub
  end function calc_magne_sum_ising2d_gpu
end module ising2d_GPU_m
