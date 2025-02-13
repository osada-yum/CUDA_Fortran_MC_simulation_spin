program xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs
  use, intrinsic :: iso_fortran_env
  use mpi
  use output_utilities_m
  use xy2d_periodic_gpu_m
  use variance_kahan_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: outs(*) = [output_unit, error_unit]
  integer(int32), parameter :: mcs = 100000_int32
  integer(int32), parameter :: tot_sample = 2222_int32
  integer(int64), parameter :: nx = 1500_int64
  integer(int64), parameter :: ny = nx
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 0.890d0
  integer(int32), parameter :: iseed = 42_int32
  integer(int64), parameter :: n_skip = 0_int64
  type(xy2d_gpu) :: xy2d
  real(real64) :: mx, my, mabs, e, ac
  type(variance_covariance_kahan) :: order_parameter_abs(mcs)
  type(variance_kahan) :: order_parameter_x(mcs), order_parameter_y(mcs), autocorrelation(mcs)
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
     write(outs(i), '(*(g0, 1x))' ) "# Fix the magnetization at 1MCS as x-aligned"
  end do

  do sample = 1, tot_sample
     call xy2d%set_random_spin()
     call xy2d%set_initial_magne_autocorrelation_state()

     write(error_unit, '(*(g0, 1x))') "#", sample, xy2d%calc_magne_sum() * n_inv_r64
     write(error_unit, '(*(a, i0))') "# Sample: ", sample, " / ", tot_sample

     do i = 1, mcs
        call xy2d%update()
        if (i == 1) call xy2d%rotate_summation_magne_and_autocorrelation_toward_xaxis()
        mx = xy2d%calc_magne_sum()
        my = xy2d%calc_magne_y_sum()
        call order_parameter_x(i)%add_data(mx * n_inv_r64)
        call order_parameter_y(i)%add_data(my * n_inv_r64)
        e = xy2d%calc_energy_sum()
        ! write(error_unit, '(*(g0, 1x))') i, mx, e
        mabs = hypot(mx * n_inv_r64, my * n_inv_r64)
        call order_parameter_abs(i)%add_data(mabs, e * n_inv_r64)
        ac = xy2d%calc_autocorrelation_sum()
        call autocorrelation(i)%add_data(ac * n_inv_r64)
        ! write(error_unit, *) sample, i, mx * n_inv_r64
     end do
  end do

  call output_abs_parameters_from_disorder(xy2d%nall(), mcs, order_parameter_abs, order_parameter_x, order_parameter_y, autocorrelation)
end program xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs
