!> CUDAFortran implementation for 2-dimensional Ising model.
module clock_gpu_multi_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  integer(int32), public, protected :: clock_gpu_stat
  integer(int32), parameter :: NUM_THREADS_X = 512, NUM_THREADS_Y = 2
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
     integer(int32) :: n_multi_
     real(real64), allocatable, device :: spin_magne_(:)
     integer(int32), allocatable, device :: spins_(:, :)
     real(real64), allocatable, device :: energy_table_(:, :, :)
     real(real64), allocatable, device :: ws_(:, :, :, :, :, :)
     real(real64), allocatable, device :: randoms_(:, :), next_states_(:, :)
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
  impure subroutine init_clock_gpu(this, nx, ny, kbt, state, n_multi)
    class(clock_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: state, n_multi
    integer(int32) :: i
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%max_state_ = state
    this%pi_state_inv_ = 2 * pi / state
    this%norishiro_begin_ = 1 - nx
    this%norishiro_end_ = this%nall_ + nx

    this%n_multi_ = n_multi
    allocate(this%spins_(this%norishiro_begin_:this%norishiro_end_, n_multi))
    allocate(this%randoms_(1:this%nall_, n_multi))
    allocate(this%next_states_(1:this%nall_, n_multi))
    block
      real(real64) :: spin_magne(0:state-1)
      do i = 0, state - 1
         spin_magne(i) = cos(this%pi_state_inv_ * i)
      end do
      allocate(this%spin_magne_(0:state-1), source = spin_magne)
    end block
    clock_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    call this%set_allup_spin()
    allocate(this%energy_table_(0:state-1, 0:state-1, 0:state-1))
    allocate(this%ws_(0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1))
    call this%set_kbt(kbt)
  end subroutine init_clock_gpu
  !> set_allup_spin_clock_gpu: Set 'ferromagnetic' initial state of Ising.
  pure subroutine set_allup_spin_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    this%spins_(:, :) = 0_int32
  end subroutine set_allup_spin_clock_gpu
  !> set_random_spin_clock_gpu: Set 'paramagnetic' initial state of Ising.
  impure subroutine set_random_spin_clock_gpu(this)
    class(clock_gpu), intent(inout) :: this
    type(dim3) :: blocks, threads
    blocks = dim3((this%nall_ + NUM_THREADS_X - 1) / NUM_THREADS_X, (this%n_multi_ + NUM_THREADS_Y - 1) / NUM_THREADS_Y, 1)
    threads = dim3(NUM_THREADS_X, NUM_THREADS_Y, 1)
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_ * this%n_multi_)
    call set_random_spin_sub <<<blocks, threads>>>&
         & (this%nall_, this%n_multi_, this%spins_(1:this%nall_, 1:this%n_multi_), this%randoms_(1:this%nall_, 1:this%n_multi_), this%max_state_)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine set_random_spin_clock_gpu
  attributes(global) impure subroutine set_random_spin_sub(n, n_multi, spins, randoms, max_state)
    integer(int64), value :: n
    integer(int32), value :: n_multi
    integer(int32), intent(inout) :: spins(n, n_multi)
    real(real64), intent(in) :: randoms(n, n_multi)
    integer(int32), value :: max_state
    integer(int64) :: idx_x, idx_y
    idx_x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    idx_y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (idx_x > n .or. idx_y > n_multi) return !> over norishiro.
    !> 1 <= idx <= this%nall_.
    spins(idx_x, idx_y) = floor(randoms(idx_x, idx_y) * max_state)
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
    type(dim3) :: blocks, threads
    blocks = dim3((this%nall_ + NUM_THREADS_X - 1) / NUM_THREADS_X, (this%n_multi_ + NUM_THREADS_Y - 1) / NUM_THREADS_Y, 1)
    threads = dim3(NUM_THREADS_X, NUM_THREADS_Y, 1)
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    call update_norishiro_sub <<<blocks, threads>>> &
         & (lb, ub, this%nx_, this%nall_, this%n_multi_, this%spins_)
    clock_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_clock_gpu
  attributes(global) pure subroutine update_norishiro_sub(lb, ub, nx, nall, n_multi, spins)
    integer(int64), value :: lb, ub, nx, nall
    integer(int32), value :: n_multi
    integer(int32), intent(inout) :: spins(lb:ub, n_multi)
    integer(int64) :: idx_x, idx_y
    idx_x = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    idx_y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (idx_x > nx .or. idx_y > n_multi) return !> over norishiro.
    !> 1 <= idx <= this%nx_
    ! norishiro top, [nall+1:nall+nx] <- [1:nx]
    spins(nall + idx_x, idx_y) = spins(idx_x, idx_y)
    ! norishiro bottom, [-nx:0] <- [nall-nx+1:nall]
    spins(idx_x - nx, idx_y) = spins(nall - nx + idx_x, idx_y)
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
    type(dim3) :: blocks, threads
    blocks = dim3((this%nall_ + NUM_THREADS_X - 1) / NUM_THREADS_X, (this%n_multi_ + NUM_THREADS_Y - 1) / NUM_THREADS_Y, 1)
    threads = dim3(NUM_THREADS_X, NUM_THREADS_Y, 1)
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_ * this%n_multi_)
    clock_gpu_stat = curandGenerate(this%rand_gen_, this%next_states_, this%nall_ * this%n_multi_)
    call update_sub <<<blocks, threads>>> &
         & (lb, ub, this%nx_, this%nall_, this%n_multi_, this%spins_, this%randoms_, this%next_states_, this%max_state_, this%ws_, 1)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    call update_sub <<<blocks, threads>>> &
         & (lb, ub, this%nx_, this%nall_, this%n_multi_, this%spins_, this%randoms_, this%next_states_, this%max_state_, this%ws_, 2)
    clock_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine update_clock_gpu
  attributes(global) pure subroutine update_sub(lb, ub, nx, nall, n_multi, spins, randoms, next_states, max_state, ws, offset)
    integer(int64), value :: lb, ub, nx, nall
    integer(int32), value :: n_multi
    integer(int32), intent(inout) :: spins(lb:ub, n_multi)
    real(real64), intent(in) :: randoms(nall, n_multi), next_states(nall, n_multi)
    real(real64), intent(in) :: ws(0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1, 0:max_state-1)
    integer(int32), value :: max_state
    integer(int32) :: next_state
    integer(int32), value :: offset
    integer(int64) :: idx_x, idx_y
    idx_x = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 2 + offset
    idx_y = (blockIdx%y - 1) * blockDim%y + threadIdx%y
    if (idx_x > nall .or. idx_y > n_multi) return
    !> do idx_x = offset, this%nall_, 2
    !> do idx_y = 1, this%n_multi_
    next_state = floor(next_states(idx_x, idx_y) * max_state)
    if (randoms(idx_x, idx_y) >= ws(spins(idx_x + nx, idx_y), spins(idx_x - nx, idx_y), &
         & spins(idx_x - 1, idx_y), spins(idx_x + 1, idx_y), spins(idx_x, idx_y), next_state)) &
         & return
    !> randoms(idx_x, idx_y) < exparr(delta_energy)
    spins(idx_x, idx_y) = next_state
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
    integer(int32), allocatable :: res(:, :)
    allocate(res, source = this%spins_)
  end function spins_clock_gpu

  !> calc_energy_sum_clock_gpu: Calculate summation of energy.
  pure subroutine calc_energy_sum_clock_gpu(this, res)
    class(clock_gpu), intent(in) :: this
    real(real64), intent(inout) :: res(:)
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    res(:) = calc_energy_sum_sub(lb, ub, this%nall_, this%max_state_, this%n_multi_, this%spins_(:, :), this%energy_table_)
  contains
    pure function calc_energy_sum_sub(lb, ub, n, max_state, n_multi, spins, energy_table) result(res)
      integer(int64), intent(in) :: lb, ub, n
      integer(int32), intent(in) :: max_state, n_multi
      integer(int32), intent(in), device :: spins(lb:ub, n_multi)
      real(real64), intent(in), device :: energy_table(0:max_state-1, 0:max_state-1, 0:max_state-1)
      real(real64) :: res(n_multi)
      integer(int64) :: i
      integer(int32) :: j
      res(:) = 0.0_real64
      do j = 1, n_multi
         !$acc parallel loop present(spins, energy_table) reduction(+:res)
         do i = 1, n
            res(j) = res(j) + energy_table(spins(i - this%nx_, j), spins(i - 1, j), spins(i, j))
         end do
      end do
    end function calc_energy_sum_sub
  end subroutine calc_energy_sum_clock_gpu
  !> calc_magne_sum_clock_gpu: Calculate summation of energy.
  pure subroutine calc_magne_sum_clock_gpu(this, res)
    class(clock_gpu), intent(in) :: this
    real(real64), intent(inout) :: res(:)
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    res(:) = calc_magne_sum_sub(lb, ub, this%nall_, this%max_state_, this%n_multi_, this%spins_(:, :), this%spin_magne_)
  contains
    pure function calc_magne_sum_sub(lb, ub, n, max_state, n_multi, spins, spin_magne) result(res)
      integer(int64), intent(in) :: lb, ub, n
      integer(int32), intent(in) :: max_state, n_multi
      integer(int32), intent(in), device :: spins(lb, ub, n_multi)
      real(real64), intent(in), device :: spin_magne(0:max_state-1)
      real(real64) :: res(n_multi)
      integer(int64) :: i
      integer(int32) :: j
      res(:) = 0.0_real64
      do j = 1, n_multi
      !$acc parallel loop present(spins, spin_magne) reduction(+:res)
         do i = 1, n
            res(j) = res(j) + spin_magne(spins(i, j))
         end do
      end do
    end function calc_magne_sum_sub
  end subroutine calc_magne_sum_clock_gpu
end module clock_gpu_multi_m
