program clock_simple_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use clock_simple_gpu_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: outs(*) = [output_unit, error_unit]

  integer(int32), parameter :: mcs = 10000_int32
  integer(int32), parameter :: tot_sample = 100_int32

  integer(int32), parameter :: iseed = 42_int32
  integer(int64), parameter :: n_skip = 0_int64

  real(real64) :: m, e
  type(variance_covariance_kahan), allocatable :: order_parameter(:)
  integer(int32) :: i, sample

  call init_sixclock(iseed)
  call skip_curand_clock(2 * n_skip * nx * ny * (mcs + 1) * tot_sample)

  do i = 1, size(outs)
     write(outs(i), '(a, i0)') "# size: ", nall
     write(outs(i), '(a, 2(i0, 1x))') "# nx, ny: ", nx, ny
     write(outs(i), '(a, i0)') "# state: ", mstate
     write(outs(i), '(a, i0)') "# sample: ", tot_sample
     write(outs(i), '(a, i0)') "# mcs: ", mcs
     write(outs(i), '(a, g0)') "# kbt: ", kbt
     write(outs(i), '(a, i0)' ) "# initial seed: ", iseed
     write(outs(i), '(a)' ) "# method: Metropolis(GPU)"
  end do

  allocate(order_parameter(mcs), source = variance_covariance_kahan())
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call init_sixclock_order()
     do i = 1, mcs
        call update_metropolis()
        m = calc_magne()
        e = calc_energy()
        ! e = clock%calc_energy_sum()
        call order_parameter(i)%add_data(m, e)
     end do
  end do
  write(output_unit, '(g0)') "# Nsize Nsample time <m> <e> <me> <m^2> <e^2> χ C m'"
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') nall, order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), order_parameter(i)%mean_v1v2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & nall * order_parameter(i)%var1(), nall * order_parameter(i)%var2(), &
          & nall * order_parameter(i)%cov()
  end do
end program clock_simple_gpu_relaxation
