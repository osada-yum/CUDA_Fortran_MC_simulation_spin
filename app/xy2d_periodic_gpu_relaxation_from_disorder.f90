program xy2d_periodic_gpu_relaxation_from_disorder
  use, intrinsic :: iso_fortran_env
  use xy2d_periodic_gpu_m
  use variance_kahan_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: outs(*) = [output_unit, error_unit]
  integer(int32), parameter :: mcs = 100_int32
  integer(int32), parameter :: tot_sample = 50000_int32
  integer(int64), parameter :: nx = 1000_int64
  integer(int64), parameter :: ny = nx
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 0.890d0
  integer(int32), parameter :: iseed = 42_int32
  integer(int64), parameter :: n_skip = 0_int64
  type(xy2d_gpu) :: xy2d
  real(real64) :: m, e, ac
  type(variance_covariance_kahan) :: order_parameter(mcs)
  type(variance_kahan) :: autocorrelation(mcs)
  integer(int32) :: i, sample

  call xy2d%init(nx, ny, kbt, iseed)
  call xy2d%skip_curand(2 * n_skip * nx * ny * (mcs + 1) * tot_sample)

  do i = 1, size(outs)
     write(outs(i), '(a, i0)') "# size: ", xy2d%nall()
     write(outs(i), '(a, i0, 1x, i0)') "# nx, ny: ", xy2d%nx(), xy2d%ny()
     write(outs(i), '(a, i0)') "# sample: ", tot_sample
     write(outs(i), '(a, i0)') "# mcs: ", mcs
     write(outs(i), '(a, g0)') "# kbt: ", kbt
     write(outs(i), '(a, i0)' ) "# initial seed: ", iseed
     write(outs(i), '(a, i0)' ) "# n_skip seed: ", n_skip
     write(outs(i), '(a)' ) "# method: Metropolis"
     write(outs(i), '(a)' ) "# initial state: disorder"
  end do

  do sample = 1, tot_sample
     call xy2d%set_random_spin()
     call xy2d%set_initial_magne_autocorrelation_state()

     write(error_unit, '(*(g0, 1x))') "#", sample, xy2d%calc_magne_sum() * n_inv_r64
     write(error_unit, '(*(a, i0))') "# Sample: ", sample, " / ", tot_sample

     do i = 1, mcs
        call xy2d%update()
        m = xy2d%calc_magne_sum()
        e = xy2d%calc_energy_sum()
        ! write(error_unit, '(*(g0, 1x))') i, m, e
        call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
        ac = xy2d%calc_autocorrelation_sum()
        call autocorrelation(i)%add_data(ac * n_inv_r64)
        ! write(error_unit, *) sample, i, m * n_inv_r64
     end do
  end do
  write(output_unit, '(*(g0, 1x))') "# N, Nsample, time, <m>, <e>, <m^2>, <e^2>, N*Var[m], N*Var[e], N*Cov[m,e], <A>, <A^2>, N*Var[A]"
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') xy2d%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & xy2d%nall() * order_parameter(i)%var1(), xy2d%nall() * order_parameter(i)%var2(), &
          & xy2d%nall() * order_parameter(i)%cov(), &
          & autocorrelation(i)%mean(), &
          & autocorrelation(i)%square_mean(), &
          & xy2d%nall() * autocorrelation(i)%var()
  end do
end program xy2d_periodic_gpu_relaxation_from_disorder
