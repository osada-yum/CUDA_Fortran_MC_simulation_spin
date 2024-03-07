program clock_gpu_multi_relaxation
  use, intrinsic :: iso_fortran_env
  use clock_gpu_multi_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 100000_int32
  integer(int32), parameter :: tot_sample = 150_int32
  integer(int64), parameter :: nx = 501_int64
  integer(int64), parameter :: ny = 500_int64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 0.80_real64
  integer(int32), parameter :: state = 6_int32
  integer(int32), parameter :: n_multi = 2_int32
  type(clock_gpu) :: clock
  real(real64) :: m(n_multi), e(n_multi)
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, j, sample
  call clock%init(nx, ny, kbt, state, n_multi)

  write(error_unit, '(a, i0)') "# size: ", clock%nall()
  write(error_unit, '(a, 3(i0, 1x))') "# nx, ny: ", clock%nx(), clock%ny(), state
  write(error_unit, '(2(a, i0))') "# sample: ", tot_sample, " x ", n_multi
  write(error_unit, '(a, i0)') "# mcs: ", mcs
  write(error_unit, '(a, g0)') "# kbt: ", kbt
  write(error_unit, '(a)' ) "# method: Metropolis"

  write(output_unit, '(a, i0)') "# size: ", clock%nall()
  write(output_unit, '(a, 3(i0, 1x))') "# nx, ny: ", clock%nx(), clock%ny(), state
  write(output_unit, '(2(a, i0))') "# sample: ", tot_sample, " x ", n_multi
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  write(output_unit, '(a)' ) "# method: Metropolis"
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call clock%set_allup_spin()
     ! call clock%set_random_spin()
     do i = 1, mcs
        call clock%update()
        call clock%calc_magne_sum(m)
        call clock%calc_energy_sum(e)
        ! write(error_unit, '(*(g0, 1x))') i, m, e
        do j = 1, n_multi
           call order_parameter(i)%add_data(m(j) * n_inv_r64, e(j) * n_inv_r64)
        end do
     end do
  end do
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') clock%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & order_parameter(i)%var1(), order_parameter(i)%var2(), &
          & order_parameter(i)%cov()
  end do
end program clock_gpu_multi_relaxation
