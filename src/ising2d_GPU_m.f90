module ising2d_GPU_m
  use, intrinsic :: iso_fortran_env
  implicit none
  type :: ising2d_GPU
     integer(int64) :: nx_, ny_
     integer(int8), allocatable :: ising_(:)
  end type ising2d_GPU
contains
end module ising2d_GPU_m
