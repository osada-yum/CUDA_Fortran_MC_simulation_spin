program ising2d_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use ising2d_gpu_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000_int32
  integer(int32), parameter :: tot_sample = 1440000_int32
  integer(int64), parameter :: nx = 1001_int64
  integer(int64), parameter :: ny = 1000_int64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 2.26918531421_real64
  integer(int32), parameter :: iseed = 42
  type(ising2d_gpu) :: ising2d
  integer(int64) :: m, e
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, sample
  call ising2d%init(nx, ny, kbt, iseed)
  write(error_unit, '(a, i0)') "# size: ", ising2d%nall()
  write(error_unit, '(a, i0, 1x, i0)') "# nx, ny: ", ising2d%nx(), ising2d%ny()
  write(error_unit, '(a, i0)') "# sample: ", tot_sample
  write(error_unit, '(a, i0)') "# mcs: ", mcs
  write(error_unit, '(a, g0)') "# kbt: ", kbt
  write(error_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(error_unit, '(a)' ) "# method: Metropolis"

  write(output_unit, '(a, i0)') "# size: ", ising2d%nall()
  write(output_unit, '(a, i0, 1x, i0)') "# nx, ny: ", ising2d%nx(), ising2d%ny()
  write(output_unit, '(a, i0)') "# sample: ", tot_sample
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  write(output_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(output_unit, '(a)' ) "# method: Metropolis"

  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call ising2d%set_allup_spin()
     ! call ising2d%set_random_spin()
     do i = 1, mcs
        call ising2d%update()
        m = ising2d%calc_magne_sum()
        e = ising2d%calc_energy_sum()
        call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
     end do
  end do
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') ising2d%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & ising2d%nall() * order_parameter(i)%var1(), ising2d%nall() * order_parameter(i)%var2(), &
          & ising2d%nall() * order_parameter(i)%cov()
  end do
end program ising2d_gpu_relaxation
