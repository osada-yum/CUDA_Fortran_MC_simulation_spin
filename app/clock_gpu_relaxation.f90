program clock_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use clock_gpu_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 100000_int32
  integer(int32), parameter :: tot_sample = 100_int32
  integer(int64), parameter :: nx = 501_int64
  integer(int64), parameter :: ny = 500_int64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 0.80_real64
  integer(int32), parameter :: state = 6_int32
  integer(int32), parameter :: iseed = 42
  type(clock_gpu) :: clock
  real(real64) :: m, e
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, sample
  call clock%init(nx, ny, kbt, state, iseed)
  write(error_unit, '(a, i0)') "# size: ", clock%nall()
  write(error_unit, '(a, 3(i0, 1x))') "# nx, ny: ", clock%nx(), clock%ny(), state
  write(error_unit, '(a, i0)') "# sample: ", tot_sample
  write(error_unit, '(a, i0)') "# mcs: ", mcs
  write(error_unit, '(a, g0)') "# kbt: ", kbt
  write(error_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(error_unit, '(a)' ) "# method: Metropolis"

  write(output_unit, '(a, i0)') "# size: ", clock%nall()
  write(output_unit, '(a, 3(i0, 1x))') "# nx, ny: ", clock%nx(), clock%ny(), state
  write(output_unit, '(a, i0)') "# sample: ", tot_sample
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  write(output_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(output_unit, '(a)' ) "# method: Metropolis"
  allocate(order_parameter(mcs), source = variance_covariance_kahan())
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call clock%set_allup_spin()
     ! call clock%set_random_spin()
     do i = 1, mcs
        call clock%update()
        m = clock%calc_magne_sum()
        e = 0.0_real64
        ! e = clock%calc_energy_sum()
        call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
     end do
  end do
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') clock%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & clock%nall() * order_parameter(i)%var1(), clock%nall() * order_parameter(i)%var2(), &
          & clock%nall() * order_parameter(i)%cov()
  end do
end program clock_gpu_relaxation
