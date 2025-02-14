program xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs
  use, intrinsic :: iso_fortran_env
  use xy2d_periodic_gpu_m
  implicit none
  integer(int64), parameter :: nx = 2_int64**2
  integer(int64), parameter :: ny = nx
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  real(real64), parameter :: kbt = 0.890d0
  integer(int32), parameter :: iseed = 42_int32
  type(xy2d_gpu) :: xy2d

  call xy2d%init(nx, ny, kbt, iseed)

  call xy2d%set_allup_spin()

  write(error_unit, '(*(g0, 1x))') "#", xy2d%calc_energy_sum() * n_inv_r64
  write(error_unit, '(*(g0, 1x))') "#", xy2d%calc_magne_sum() * n_inv_r64
end program xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs
