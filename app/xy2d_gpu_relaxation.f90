program xy2d_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use xy2d_gpu_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000_int32
  integer(int32), parameter :: tot_sample = 1440000_int32
  integer(int64), parameter :: nx = 1001_int64
  integer(int64), parameter :: ny = 1000_int64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 1.0_real64
  integer(int32), parameter :: iseed = 42
  type(xy2d_gpu) :: xy2d
  real(real64) :: m, e
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, sample
  call xy2d%init(nx, ny, kbt, iseed)
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call xy2d%set_allup_spin()
     ! call xy2d%set_random_spin()
     do i = 1, mcs
        call xy2d%update()
        m = xy2d%calc_magne_sum()
        e = xy2d%calc_energy_sum()
        call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
     end do
  end do
  write(output_unit, '(a, i0)') "# size: ", xy2d%nall()
  write(output_unit, '(a, i0, 1x, i0)') "# nx, ny: ", xy2d%nx(), xy2d%ny()
  write(output_unit, '(a, i0)') "# sample: ", tot_sample
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  write(output_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(output_unit, '(a)' ) "# method: Metropolis"
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') xy2d%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & order_parameter(i)%var1(), order_parameter(i)%var2(), &
          & order_parameter(i)%cov()
  end do
end program xy2d_gpu_relaxation
