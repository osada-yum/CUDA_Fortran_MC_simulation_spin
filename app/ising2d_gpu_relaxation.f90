program ising2d_gpu_relaxation
  use, intrinsic :: iso_fortran_env
  use ising2d_gpu_m
  implicit none
  integer(int32), parameter :: mcs = 1000_int32, tot_sample = 100_int32
  integer(int64), parameter :: nx = 1001_int64, ny = 1000_int64
  real(real64), parameter :: kbt = 2.2_real64
  type(ising2d_gpu) :: ising2d
  integer(int64) :: tot_magne(mcs), tot_energy(mcs)
  integer(int32) :: i, sample
  call ising2d%init(nx, ny, kbt)
  tot_magne(:) = 0_int64
  tot_energy(:) = 0_int64
  do sample = 1, tot_sample
     write(error_unit, '(*(a, i0))') "Sample: ", sample, " / ", tot_sample
     call ising2d%set_allup_spin()
     ! call ising2d%set_random_spin()
     do i = 1, mcs
        call ising2d%update()
        tot_magne(i) = tot_magne(i) + ising2d%calc_magne_sum()
        tot_energy(i) = tot_energy(i) + ising2d%calc_energy_sum()
     end do
  end do
  write(output_unit, '(a, i0)') "# size: ", ising2d%nall()
  write(output_unit, '(a, i0)') "# sample: ", tot_sample
  write(output_unit, '(a, i0)') "# mcs: ", mcs
  write(output_unit, '(a, g0)') "# kbt: ", kbt
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') i, &
          & real(tot_magne(i), real64) / (ising2d%nall() * tot_sample), &
          & real(tot_energy(i), real64) / (ising2d%nall() * tot_sample)
  end do
end program ising2d_gpu_relaxation
