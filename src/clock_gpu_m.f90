!> CUDAFortran implementation for 2-dimensional Ising model.
module clock_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  integer(int32), public, protected :: clock_gpu_stat
  integer(int64), parameter :: NUM_THREADS = 512
  integer(int32), parameter :: clock_max_state_limit = 50 ! 実装的に 50^5 が限界.
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  public :: clock_gpu
  type :: clock_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nall_
     integer(int64) :: norishiro_begin_, norishiro_end_
     integer(int32) :: max_state_
     real(real64) :: pi_state_inv_
     real(real64), allocatable, device :: spin_magne_(:)
     integer(int32), allocatable, device :: spins_(:)
     real(real64), allocatable, device :: energy_table_(:, :, :)
     real(real64), allocatable, device :: ws_(:, :, :, :, :, :)
     real(real64), allocatable, device :: randoms_(:), next_states_(:)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_clock_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_clock_gpu
     procedure, pass :: set_random_spin => set_random_spin_clock_gpu
     procedure, pass :: set_kbt => set_kbt_clock_gpu
     procedure, pass :: set_beta => set_beta_clock_gpu
     !> updater.
     procedure, pass, private :: update_norishiro => update_norishiro_clock_gpu
     procedure, pass, private :: update_ws => update_ws_clock_gpu
     procedure, pass :: update => update_clock_gpu
     !> getter.
     procedure, pass :: nx => nx_clock_gpu
     procedure, pass :: ny => ny_clock_gpu
     procedure, pass :: nall => nall_clock_gpu
     procedure, pass :: kbt => kbt_clock_gpu
     procedure, pass :: beta => beta_clock_gpu
     procedure, pass :: spins => spins_clock_gpu
     !> calculator.
     procedure, pass :: calc_energy_sum => calc_energy_sum_clock_gpu
     procedure, pass :: calc_magne_sum => calc_magne_sum_clock_gpu
  end type clock_gpu
contains
  impure subroutine init_clock_gpu(this, nx, ny, kbt, state, iseed)
    class(clock_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: state
    integer(int32), intent(in) :: iseed
    integer(int32) :: i
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%max_state_ = state
    this%pi_state_inv_ = 2 * pi / state
    this%norishiro_begin_ = 1 - nx
    this%norishiro_end_ = this%nall_ + nx
    allocate(this%spins_(this%norishiro_begin_:this%norishiro_end_))
    allocate(this%randoms_(1:this%nall_))
    allocate(this%next_states_(1:this%nall_))
    block
      real(real64) :: spin_magne(0:state-1)
      do i = 0, state - 1
         spin_magne(i) = cos(this%pi_state_inv_ * i)
      end do
      allocate(this%spin_magne_(0:state-1), source = spin_magne)
    end block
    clock_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    clock_gpu_stat = curandSetPseudoRandomGeneratorSeed(this%rand_gen_, iseed)
    call this%set_allup_spin()
    allocate(this%energy_table_(0:state-1, 0:state-1, 0:state-1))
    allocate(this%ws_(0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1))
    call this%set_kbt(kbt)
  end subroutine init_clock_gpu
  !> set_allup_spin_clock_gpu: Set 'ferromagnetic' initial state of Ising.
  pure subroutine set_allup_spin_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    this%spins_(:) = 0_int32
  end subroutine set_allup_spin_clock_gpu
  !> set_random_spin_clock_gpu: Set 'paramagnetic' initial state of Ising.
  impure subroutine set_random_spin_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub <<<(this%nall_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>>&
         & (this%nall_, this%spins_(1:this%nall_), this%randoms_(1:this%nall_), this%max_state_)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine set_random_spin_clock_gpu
  attributes(global) impure subroutine set_random_spin_sub(n, spins, randoms, max_state)
    integer(int64), value :: n
    integer(int32), intent(inout) :: spins(n)
    real(real64), intent(in) :: randoms(n)
    integer(int32), value :: max_state
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > n) return !> over norishiro.
    !> 1 <= idx <= this%nall_.
    spins(idx) = floor(randoms(idx) * max_state)
  end subroutine set_random_spin_sub
  pure subroutine update_ws_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    integer(int32) :: center, center_before, center_after, i, j, k, l
    real(real64), allocatable :: energy_table(:, :, :)
    real(real64), allocatable :: ws(:, :, :, :, :, :)
    real(real64) :: delta_e
    allocate(energy_table(0:this%max_state_-1, 0:this%max_state_-1, 0:this%max_state_-1))
    do center = 0, this%max_state_ - 1
       do j = 0, this%max_state_ - 1
          do i = 0, this%max_state_ - 1
             energy_table(i, j, center) = calc_local_energy(i, j, center)
          end do
       end do
    end do
    this%energy_table_(:, :, :) = energy_table(:, :, :)
    allocate(ws(0:this%max_state_-1, 0:this%max_state_-1, 0:this%max_state_-1, 0:this%max_state_-1, 0:this%max_state_-1, 0:this%max_state_-1))
    do center_after = 0, this%max_state_ - 1
       do center_before = 0, this%max_state_ - 1
          do l = 0, this%max_state_ - 1
             do k = 0, this%max_state_ - 1
                do j = 0, this%max_state_ - 1
                   do i = 0, this%max_state_ - 1
                      delta_e = (energy_table(i, j, center_after) + energy_table(k, l, center_after)) &
                           & - (energy_table(i, j, center_before) + energy_table(k, l, center_before))
                      if (delta_e <= 0.0_real64) then
                         ws(i, j, k, l, center_before, center_after) = 1.0_real64
                      else
                         ws(i, j, k, l, center_before, center_after) = exp(- this%beta_ * delta_e)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    this%ws_(:, :, :, :, :, :) = ws(:, :, :, :, :, :)
  contains
    pure real(real64) function calc_local_energy(i, j, center) result(res)
      integer(int32), intent(in) :: i, j, center
      res = - sum(cos(this%pi_state_inv_ * [i - center, j - center]))
    end function calc_local_energy
  end subroutine update_ws_clock_gpu
  !> update_norishiro_clock_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    call update_norishiro_sub <<<(this%nx_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_)
    clock_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_clock_gpu
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
  !> set_kbt_clock_gpu: Set new kbt.
  pure subroutine set_kbt_clock_gpu(this, kbt)
    class(clock_gpu), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_clock_gpu
  !> set_beta_clock_gpu: Set new kbt.
  pure subroutine set_beta_clock_gpu(this, beta)
    class(clock_gpu), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
    call this%update_ws()
  end subroutine set_beta_clock_gpu
  !> update_clock_gpu: Update by Metropolis method.
  impure subroutine update_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_(1:this%nall_), this%nall_)
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%next_states_(1:this%nall_), this%nall_)
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%randoms_(1:this%nall_), this%next_states_(1:this%nall_), this%max_state_, this%ws_, 1)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_, this%randoms_(1:this%nall_), this%next_states_(1:this%nall_), this%max_state_, this%ws_, 2)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine update_clock_gpu
  attributes(global) pure subroutine update_sub(lb, ub, nx, nall, spins, randoms, next_states, max_state, ws, offset)
    integer(int64), value :: lb, ub, nx, nall
    integer(int32), intent(inout) :: spins(lb:ub)
    real(real64), intent(in) :: randoms(nall), next_states(nall)
    real(real64), intent(in) :: ws(0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1)
    integer(int32), value :: max_state
    integer(int32) :: next_state
    integer(int32), value :: offset
    integer(int64) :: idx
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 2 + offset
    if (idx > nall) return
    !> do idx = offset, this%nall_, 2
    next_state = floor(next_states(idx) * max_state)
    if (randoms(idx) > ws(spins(idx + nx), spins(idx - nx), spins(idx - 1), spins(idx + 1), spins(idx), next_state)) &
         & return
    !> randoms(idx) <= exparr(delta_energy)
    spins(idx) = next_state
  end subroutine update_sub

  pure integer(int64) function nx_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = this%nx_
  end function nx_clock_gpu
  pure integer(int64) function ny_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = this%ny_
  end function ny_clock_gpu
  pure integer(int64) function nall_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = this%nall_
  end function nall_clock_gpu
  pure real(real64) function kbt_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_clock_gpu
  pure real(real64) function beta_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = this%beta_
  end function beta_clock_gpu
  pure function spins_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    integer(int32), allocatable :: res(:)
    allocate(res, source = this%spins_)
  end function spins_clock_gpu

  !> calc_energy_sum_clock_gpu: Calculate summation of energy.
  pure real(real64) function calc_energy_sum_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = calc_energy_sum_sub(this%nall_, this%max_state_, this%norishiro_end_, this%spins_(1:this%norishiro_end_), this%energy_table_)
  contains
    pure real(real64) function calc_energy_sum_sub(n, max_state, norishiro_end, spins, energy_table) result(res)
      integer(int64), intent(in) :: n, norishiro_end
      integer(int32), intent(in) :: max_state
      integer(int32), intent(in), device :: spins(norishiro_end)
      real(real64), intent(in), device :: energy_table(0:max_state-1, 0:max_state-1, 0:max_state-1)
      integer(int64) :: i
      res = 0.0_real64
      !$acc parallel loop present(spins, energy_table) reduction(+:res)
      do i = 1, n
         res = res + energy_table(spins(i - this%nx_), spins(i - 1), spins(i))
      end do
      res = res
    end function calc_energy_sum_sub
  end function calc_energy_sum_clock_gpu
  !> calc_magne_sum_clock_gpu: Calculate summation of energy.
  pure real(real64) function calc_magne_sum_clock_gpu(this) result(res)
    class(clock_gpu), intent(in) :: this
    res = calc_magne_sum_sub(this%nall_, this%max_state_, this%spins_(1:this%nall_), this%spin_magne_)
  contains
    pure real(real64) function calc_magne_sum_sub(n, max_state, spins, spin_magne) result(res)
      integer(int64), intent(in) :: n
      integer(int32), intent(in) :: max_state
      integer(int32), intent(in), device :: spins(n)
      real(real64), intent(in), device :: spin_magne(0:max_state-1)
      integer(int64) :: i
      res = 0.0_real64
      !$acc parallel loop present(spins, spin_magne) reduction(+:res)
      do i = 1, n
         res = res + spin_magne(spins(i))
      end do
    end function calc_magne_sum_sub
  end function calc_magne_sum_clock_gpu
end module clock_GPU_m
