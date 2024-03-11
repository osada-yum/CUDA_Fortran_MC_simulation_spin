!> CUDAFortran implementation for 3-dimensional Ising model.
module ising3d_gpu_m
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand
  implicit none
  private
  integer(int32), public, protected :: ising3d_gpu_stat
  integer(int64), parameter :: NUM_THREADS = 512
  integer(int64), constant :: energy_table(0:3, 0:1)
  real(real64), constant :: ws(0:6, 0:1)
  integer(int32), parameter :: spin_map(0:1) = [-1, 1]
  integer(int32), constant :: spin_map_c(0:1)
  public :: ising3d_gpu
  type :: ising3d_gpu
     private
     real(real64) :: beta_
     integer(int64) :: nx_, ny_, nz_, nxy_, nall_
     integer(int64) :: norishiro_begin_, norishiro_end_
     integer(int32), allocatable, device :: spins_(:)
     real(real64), allocatable, device :: randoms_(:)
     integer(int64), allocatable :: energy_table_(:, :)
     real(real64), allocatable :: ws_(:, :)
     type(curandGenerator) :: rand_gen_
   contains
     procedure, pass :: init => init_ising3d_gpu
     procedure, pass :: skip_curand => skip_curand_ising3d_gpu
     !> setter.
     procedure, pass :: set_allup_spin => set_allup_spin_ising3d_gpu
     procedure, pass :: set_random_spin => set_random_spin_ising3d_gpu
     procedure, pass :: set_kbt => set_kbt_ising3d_gpu
     procedure, pass :: set_beta => set_beta_ising3d_gpu
     !> updater.
     procedure, pass, private :: update_ws => update_ws_ising3d_gpu
     procedure, pass, private :: update_norishiro => update_norishiro_ising3d_gpu
     procedure, pass :: update => update_ising3d_gpu
     !> getter.
     procedure, pass :: nx => nx_ising3d_gpu
     procedure, pass :: ny => ny_ising3d_gpu
     procedure, pass :: nz => nz_ising3d_gpu
     procedure, pass :: nall => nall_ising3d_gpu
     procedure, pass :: kbt => kbt_ising3d_gpu
     procedure, pass :: beta => beta_ising3d_gpu
     procedure, pass :: spins => spins_ising3d_gpu
     !> calculator.
     procedure, pass :: calc_energy_sum => calc_energy_sum_ising3d_gpu
     procedure, pass :: calc_magne_sum => calc_magne_sum_ising3d_gpu
  end type ising3d_gpu
contains
  impure subroutine init_ising3d_gpu(this, nx, ny, nz, kbt, iseed)
    class(ising3d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny, nz
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: iseed
    this%nx_ = nx
    this%ny_ = ny
    this%nz_ = nz
    this%nxy_ = nx * ny
    this%nall_ = nx * ny * nz
    this%norishiro_begin_ = 1 - this%nxy_
    this%norishiro_end_ = this%nall_ + this%nxy_
    allocate(this%spins_(this%norishiro_begin_:this%norishiro_end_))
    allocate(this%randoms_(1:this%nall_))
    ising3d_gpu_stat = curandCreateGenerator(this%rand_gen_, CURAND_RNG_PSEUDO_XORWOW)
    ising3d_gpu_stat = curandSetPseudoRandomGeneratorSeed(this%rand_gen_, iseed)
    call this%set_allup_spin()
    spin_map_c(:) = spin_map(:)
    allocate(this%energy_table_(0:3, 0:1))
    allocate(this%ws_(0:6, 0:1))
    call this%set_kbt(kbt)
  end subroutine init_ising3d_gpu
 impure subroutine skip_curand_ising3d_gpu(this, n_skip)
    class(ising3d_gpu), intent(inout) :: this
    integer(int64), intent(in) :: n_skip
    if (n_skip == 0_int64) return
    ising3d_gpu_stat = curandSetGeneratorOffset(this%rand_gen_, n_skip)
  end subroutine skip_curand_ising3d_gpu
  !> set_allup_spin_ising3d_gpu: Set 'ferromagnetic' initial state of Ising.
  pure subroutine set_allup_spin_ising3d_gpu(this)
    class(ising3d_gpu), intent(inout) :: this
    this%spins_(:) = 1_int32
  end subroutine set_allup_spin_ising3d_gpu
  !> set_random_spin_ising3d_gpu: Set 'paramagnetic' initial state of Ising.
  impure subroutine set_random_spin_ising3d_gpu(this)
    class(ising3d_gpu), intent(inout) :: this
    ising3d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_, this%nall_)
    call set_random_spin_sub <<<(this%nall_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>>(this%nall_, this%spins_(1:this%nall_), this%randoms_(1:this%nall_))
    ising3d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine set_random_spin_ising3d_gpu
  attributes(global) pure subroutine set_random_spin_sub(n, spins, randoms)
    integer(int64), value :: n
    integer(int32), intent(inout) :: spins(n)
    real(real64), intent(in) :: randoms(n)
    integer(int64) :: idx
    idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (idx > n) return !> over norishiro.
    !> 1 <= idx <= this%nall_.
    spins(idx) = merge(1, 0, randoms(idx) < 0.5_real64)
  end subroutine set_random_spin_sub
  !> update_norishiro_ising3d_gpu: Update norishiro by GPU.
  impure subroutine update_norishiro_ising3d_gpu(this)
    class(ising3d_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    call update_norishiro_sub <<<(this%nx_ + NUM_THREADS - 1)/NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nall_, this%spins_)
    ising3d_gpu_stat = cudaDeviceSynchronize()
  end subroutine update_norishiro_ising3d_gpu
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
  !> set_kbt_ising3d_gpu: Set new kbt.
  pure subroutine set_kbt_ising3d_gpu(this, kbt)
    class(ising3d_gpu), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_ising3d_gpu
  !> set_beta_ising3d_gpu: Set new kbt.
  pure subroutine set_beta_ising3d_gpu(this, beta)
    class(ising3d_gpu), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
    call this%update_ws()
  end subroutine set_beta_ising3d_gpu
  !> update_ws_ising3d_gpu: Update exparr_(:) by temperature.
  !> Because Ising model is discrete, we can calculate `exp(-βΔE)` in advance.
  pure subroutine update_ws_ising3d_gpu(this)
    class(ising3d_gpu), intent(inout) :: this
    real(real64) :: tmp
    integer(int64) :: e1, e2
    integer(int32) :: i1, i2, i3, i4, i5, i6, s1, s2
    do i1 = 0, 1
       do i2 = 0, 1
          do i3 = 0, 1
             s1 = i1 + i2 + i3
             this%energy_table_(s1, 0) = - spin_map(0) * sum(spin_map([i1, i2, i3]))
             this%energy_table_(s1, 1) = - spin_map(1) * sum(spin_map([i1, i2, i3]))
          end do
       end do
    end do
    energy_table(:, :) = this%energy_table_(:, :)
    do i1 = 0, 1
       do i2 = 0, 1
          do i3 = 0, 1
             s1 = i1 + i2 + i3
             do i4 = 0, 1
                do i5 = 0, 1
                   do i6 = 0, 1
                      s2 = i4 + i5 + i6
                      e1 = this%energy_table_(s1, 0) + this%energy_table_(s2, 0)
                      e2 = this%energy_table_(s1, 1) + this%energy_table_(s2, 1)
                      this%ws_(s1 + s2, 0) = min(1.0_real64, exp(-this%beta_ * (e2 - e1)))
                      this%ws_(s1 + s2, 1) = min(1.0_real64, exp(-this%beta_ * (e1 - e2)))
                   end do
                end do
             end do
          end do
       end do
    end do
    ws(:, :) = this%ws_(:, :)
  end subroutine update_ws_ising3d_gpu
  !> update_ising3d_gpu: Update by Metropolis method.
  impure subroutine update_ising3d_gpu(this)
    class(ising3d_gpu), intent(inout) :: this
    integer(int64) :: lb, ub
    lb = lbound(this%spins_, dim = 1, kind = int64)
    ub = ubound(this%spins_, dim = 1, kind = int64)
    ising3d_gpu_stat = curandGenerate(this%rand_gen_, this%randoms_(1:this%nall_), this%nall_)
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nxy_, this%nall_, this%spins_, this%randoms_(1:this%nall_), 1)
    ising3d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
    call update_sub <<<(this%nall_ + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>> &
         & (lb, ub, this%nx_, this%nxy_, this%nall_, this%spins_, this%randoms_(1:this%nall_), 2)
    ising3d_gpu_stat = cudaDeviceSynchronize()
    call this%update_norishiro()
  end subroutine update_ising3d_gpu
  attributes(global) pure subroutine update_sub(lb, ub, nx, nxy, nall, spins, randoms, offset)
    integer(int64), value :: lb, ub, nx, nxy, nall
    integer(int32), intent(inout) :: spins(lb:ub)
    real(real64), intent(in) :: randoms(nall)
    integer(int32), value :: offset
    integer(int32) :: sum_spin
    integer(int64) :: idx
    idx = 2 * ((blockIdx%x - 1) * blockDim%x + threadIdx%x) - 2 + offset
    if (idx > nall) return
    !> do idx = offset, this%nall_, 2
    sum_spin = &
         & spins(idx - 1) + spins(idx + 1) + &
         & spins(idx - nx) + spins(idx + nx) + &
         & spins(idx - nxy) + spins(idx + nxy)
    if (randoms(idx) >= ws(sum_spin, spins(idx))) return
    !> randoms(idx) < ws(delta_energy)
    spins(idx) = 1 - spins(idx)
  end subroutine update_sub

  pure integer(int64) function nx_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = this%nx_
  end function nx_ising3d_gpu
  pure integer(int64) function ny_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = this%ny_
  end function ny_ising3d_gpu
  pure integer(int64) function nz_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = this%nz_
  end function nz_ising3d_gpu
  pure integer(int64) function nall_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = this%nall_
  end function nall_ising3d_gpu
  pure real(real64) function kbt_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_ising3d_gpu
  pure real(real64) function beta_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = this%beta_
  end function beta_ising3d_gpu
  pure function spins_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    integer(int32), allocatable :: res(:)
    allocate(res, source = this%spins_)
  end function spins_ising3d_gpu

  !> calc_energy_sum_ising3d_gpu: Calculate summation of energy.
  pure integer(int64) function calc_energy_sum_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = calc_energy_sum_sub(this%nall_, this%norishiro_end_, this%spins_(1:this%norishiro_end_), energy_table)
  contains
    pure integer(int64) function calc_energy_sum_sub(n, norishiro_end, spins, energy_table) result(res)
      integer(int64), intent(in) :: n, norishiro_end
      integer(int32), intent(in), device :: spins(norishiro_end)
      integer(int64), intent(in), device :: energy_table(0:3, 0:1)
      integer(int64) :: i
      res = 0_int64
      !$acc data present(energy_table(0:3, 0:1), spins)
      !$acc parallel loop reduction(+:res)
      do i = 1, n
         ! res = res + (spins(i + 1) + spins(i + this%nx_) + spins(i + this%nxy_)) * spins(i)
         res = res + energy_table(spins(i + 1) + spins(i + this%nx_) + spins(i + this%nxy_), spins(i))
      end do
      !$acc end data
    end function calc_energy_sum_sub
  end function calc_energy_sum_ising3d_gpu
  !> calc_magne_sum_ising3d_gpu: Calculate summation of energy.
  pure integer(int64) function calc_magne_sum_ising3d_gpu(this) result(res)
    class(ising3d_gpu), intent(in) :: this
    res = calc_magne_sum_sub(this%nall_, this%spins_(1:this%nall_))
  contains
    pure integer(int64) function calc_magne_sum_sub(n, spins) result(res)
      integer(int64), intent(in) :: n
      integer(int32), intent(in), device :: spins(n)
      integer(int64) :: i
      res = 0_int64
      !$acc data present(spins)
      !$acc parallel loop reduction(+:res)
      do i = 1, n
         res = res + spins(i)
      end do
      !$acc end data
      res = 2 * res - n
    end function calc_magne_sum_sub
  end function calc_magne_sum_ising3d_gpu
end module ising3d_GPU_m
