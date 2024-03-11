program ising3d_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use ising3d_gpu_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000_int32
  integer(int32), parameter :: tot_sample = 1440000_int32
  integer(int64), parameter :: nx = 1001_int64
  integer(int64), parameter :: ny = 1000_int64
  integer(int64), parameter :: nz = 1000_int64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny * nz, real64)
  real(real64), parameter :: kbt = 4.51152_real64
  integer(int32), parameter :: iseed = 42
  type(ising3d_gpu) :: ising3d
  integer(int64) :: m, e
  type(variance_covariance_kahan), allocatable :: order_parameter(:)
  integer(int32) :: i, sample
  call ising3d%init(nx, ny, nz, kbt, iseed)
  write(error_unit, '(a, i0)') "# size: ", ising3d%nall()
  write(error_unit, '(a, 3(i0, 1x))') "# nx, ny, nz: ", ising3d%nx(), ising3d%ny(), ising3d%nz()
  write(error_unit, '(a, i0)') "# sample: ", tot_sample
  write(error_unit, '(a, i0)') "# mcs: ", mcs
  write(error_unit, '(a, g0)') "# kbt: ", kbt
  write(error_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(error_unit, '(a)' ) "# method: Metropolis"

  write(output_unit, '(a, i0)') "# size: ", ising3d%nall()
  write(output_unit, '(a, 3(i0, 1x))') "# nx, ny, nz: ", ising3d%nx(), ising3d%ny(), ising3d%nz()
  write(output_unit, '(a, i0)') "# sample: ", tot_sample
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  write(output_unit, '(a, i0)' ) "# initial seed: ", iseed
  write(output_unit, '(a)' ) "# method: Metropolis"

  allocate(order_parameter(mcs), source = variance_covariance_kahan())
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call ising3d%set_allup_spin()
     ! call ising3d%set_random_spin()
     do i = 1, mcs
        call ising3d%update()
        m = ising3d%calc_magne_sum()
        e = ising3d%calc_energy_sum()
        call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
     end do
  end do
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') ising3d%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & order_parameter(i)%var1(), order_parameter(i)%var2(), &
          & order_parameter(i)%cov()
  end do
end program ising3d_gpu_relaxation
