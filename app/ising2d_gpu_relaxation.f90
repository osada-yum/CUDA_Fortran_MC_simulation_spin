program ising2d_gpu_relaxation
  use ising2d_gpu_m
  implicit none
  type(ising2d_gpu) :: ising2d
  integer(int32), parameter :: mcs = 1000, tot_sample = 10
  integer(int64), parameter :: nx = 1001_int64, ny = 1000_int64
  real(real64), parameter :: kbt = 2.0_real64
  integer(int64) :: tot_magne(mcs), tot_energy(mcs)
  integer(int32) :: i, sample
  call ising2d%init(nx, ny, kbt)
  tot_magne(:) = 0_int64
  tot_energy(:) = 0_int64
  do sample = 1, tot_sample
     call ising2d%set_allup_spin()
     do i = 1, mcs
        call ising2d%update()
        tot_magne(i) = tot_magne(i) + ising2d%calc_magne_sum()
        tot_energy(i) = tot_energy(i) + ising2d%calc_energy_sum()
     end do
  end do
  do sample = 1, tot_sample
     write(output_unit, '(*(g0, 1x))') sample, tot_magne(sample), tot_energy(sample)
  end do
end program ising2d_gpu_relaxation
