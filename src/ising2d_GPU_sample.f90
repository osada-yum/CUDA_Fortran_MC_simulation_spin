!!! CUDAFortran implementation for 2-dimensional Ising model.
module ising2d_GPU_m
    use, intrinsic :: iso_fortran_env
    use cudafor
    use curand
    implicit none
    private
    integer, public, parameter :: ikind = int32, rkind = real64
    integer        , parameter :: NUM_THREADS = 512

    public :: IG2_property
    type :: IG2_property
        real(rkind) :: kbt, beta
        integer     :: nx, ny, offset, particles
        integer     :: begin, end, all, all_begin, all_end
    end type IG2_property

    public :: ising2d_GPU_init
    public :: ising2d_GPU_set_kbt, ising2d_GPU_arr_set_order
    public :: ising2d_GPU_update_Metropolis
    public :: ising2d_GPU_calc_energy, ising2d_GPU_calc_magne
contains

    !! ising2d_GPU_init: Set parameters and allocate arrays.
    subroutine ising2d_GPU_init(nx, ny, kbt, prop, ising_arr, rnd_arr, exp_arr)
        integer(ikind)                  , intent(in)  :: nx, ny
        real(rkind)                     , intent(in)  :: kbt
        type(IG2_property)              , intent(out) :: prop
        integer    , allocatable, device, intent(out) :: ising_arr(:)
        real(rkind), allocatable, device, intent(out) :: rnd_arr(:), exp_arr(:)
        if (mod(nx, 2) /= 1 .and. mod(ny, 2) /= 0) then
            write(error_unit, '(a)') "nx must be odd, ny must be even"
        end if
        prop%nx = nx
        prop%ny = ny
        prop%offset    = nx
        prop%particles = nx * ny
        prop%begin     = prop%offset   + 1
        prop%end       = prop%offset   + nx*ny
        prop%all       = prop%offset*2 + nx*ny
        prop%all_begin = prop%begin-prop%offset
        prop%all_end   = prop%end  +prop%offset
        if (allocated(ising_arr)) then
            deallocate(ising_arr)
        end if
        if (allocated(rnd_arr)) then
            deallocate(rnd_arr)
        end if
        if (allocated(exp_arr)) then
            deallocate(exp_arr)
        end if
        allocate(ising_arr(prop%all_begin:prop%all_end), source = 1)
        allocate(rnd_arr(prop%begin:prop%end))
        allocate(exp_arr(-8:8))

        call ising2d_GPU_set_kbt(kbt, prop, exp_arr)
    end subroutine ising2d_GPU_init
    !! ising2d_set_kbt: set new kbt and update exparr.
    subroutine ising2d_GPU_set_kbt(kbt, prop, exp_arr)
        real(rkind)        , intent(in)    :: kbt
        type(IG2_property) , intent(inout) :: prop
        real(rkind), device, intent(out)   :: exp_arr(-8:8)
        integer                            :: i
        prop%kbt  = kbt
        prop%beta = 1/kbt
        exp_arr(:) = 1.0_rkind
        do i = 1, 8
            exp_arr(i) = exp(-prop%beta*i)
        end do
    end subroutine ising2d_GPU_set_kbt
    !! ising2d_CPU_arr_set_order: initialize ising_arr.
    subroutine ising2d_GPU_arr_set_order(prop, ising_arr)
        type(IG2_property),         intent(in)  :: prop
        integer           , device, intent(out) :: ising_arr(prop%all_begin:prop%all_end)
        ising_arr(:) = 1
    end subroutine ising2d_GPU_arr_set_order
    !! norishiro: Update norishiro.
    !                 |n|n|n|n|n| -------------^
    ! F-------------> |i|i|i|i|i| <- prop%end  |
    ! |               |i|i|i|i|i|              |
    ! |               |i|i|i|i|i|              |
    ! | prop%begin -> |i|i|i|i|i| <------------v
    ! L-------------- |n|n|n|n|n|
    attributes(global) &
    subroutine norishiro(beg, end, offset, all_beg, all_end, ising_arr)
        integer, value         :: beg, end, offset, all_beg, all_end
        integer, intent(inout) :: ising_arr(all_beg:all_end)
        integer                :: i
        i = (blockIdx%x-1) * blockDim%x + threadIdx%x
        if (i > offset) return
        ! lower norishiro.
        ising_arr(beg-offset+i-1) &
            = ising_arr(end-offset+i)
        ! upper norishiro.
        ising_arr(end+i) &
            = ising_arr(beg-offset+i-1)
    end subroutine norishiro
    !! ising2d_GPU_update_Metropolis: Update IG2_ising by Metropolis method using GPU.
    subroutine ising2d_GPU_update_Metropolis(prop, ising_arr, rnd_arr, exp_arr)
        type(IG2_property)        , intent(in)    :: prop
        integer           , device, intent(inout) :: ising_arr(prop%all_begin:prop%all_end)
        real(rkind)       , device, intent(in)    :: rnd_arr(prop%begin:prop%end), exp_arr(-8:8)
        type(dim3)                                :: dimGrid, dimBlock, dimGrid_n, dimBlock_n
        dimGrid    = dim3( (prop%all/2-1)/NUM_THREADS+1, 1, 1)   ! all sites
        dimBlock   = dim3( NUM_THREADS, 1, 1)
        dimGrid_n  = dim3( (prop%offset-1)/NUM_THREADS+1, 1, 1) ! norishiro
        dimBlock_n = dim3( NUM_THREADS, 1, 1)
        call update_checkerboard_Metropolis<<<dimGrid,dimBlock>>>&
            (0, prop%nx, prop%begin, prop%end, prop%all_begin, prop%all_end&
            , ising_arr, rnd_arr, exp_arr)
        call norishiro<<<dimGrid_n, dimBlock_n>>>&
            (prop%begin, prop%end, prop%offset, prop%all_begin, prop%all_end, ising_arr)
        call update_checkerboard_Metropolis<<<dimGrid,dimBlock>>>&
            (1, prop%nx, prop%begin, prop%end, prop%all_begin, prop%all_end&
            , ising_arr, rnd_arr, exp_arr)
        call norishiro<<<dimGrid_n, dimBlock_n>>>&
            (prop%begin, prop%end, prop%offset, prop%all_begin, prop%all_end, ising_arr)
    end subroutine ising2d_GPU_update_Metropolis
    !! update_checkerboard_Metropolis: Update sublattice of IG2_ising by Metropolis method.
    attributes(global) &
    subroutine update_checkerboard_Metropolis&
        (start_index, nx, beg, end, all_beg, all_end, ising_arr, rnd_arr, exp_arr)
        integer       , value         :: start_index, nx, beg, end, all_beg, all_end
        integer       , intent(inout) :: ising_arr(all_beg:all_end)
        real(rkind)   , intent(in)    :: rnd_arr(beg:end), exp_arr(-8:8)
        integer(int32)                :: delta_energy
        integer                       :: i
        i = 2*((blockIdx%x-1) * blockDim%x + threadIdx%x) - 2 + start_index
        if (i < beg .or. i > end) return
        delta_energy = calc_delta_energy(i, nx, all_beg, all_end, ising_arr)
        if (rnd_arr(i) < exp_arr(delta_energy)) then
            ising_arr(i) = -ising_arr(i)
        end if
    end subroutine update_checkerboard_Metropolis
    !! calc_delta_energy: Calculate delta energy of spin flip.
    ! |*|+|*|    |*|+|*|
    ! |+|+|+| -> |+|-|+|
    ! |*|+|*|    |*|+|*|
    !   -4    ->   +4   , dE = +8
    ! |*|+|*|    |*|+|*|
    ! |+|-|+| -> |+|+|+|
    ! |*|+|*|    |*|+|*|
    !   +4    ->   -4   , dE = -8
    attributes(device) &
    pure integer(int32) function calc_delta_energy(i, nx, all_beg, all_end, ising_arr) result(e)
        integer, intent(in) :: i
        integer, value      :: nx, all_beg, all_end
        integer, intent(in) :: ising_arr(all_beg:all_end)
        e = 2*ising_arr(i)*&
            ( ising_arr(i+1)&
            + ising_arr(i-1)&
            + ising_arr(i+nx)&
            + ising_arr(i-nx) )
        return
    end function calc_delta_energy
    !! ising2d_GPU_calc_energy: Calculate energy per site.
    pure real(rkind) function ising2d_GPU_calc_energy(prop, ising_arr) result(e)
        type(IG2_property)        , intent(in) :: prop
        integer           , device, intent(in) :: ising_arr(prop%all_begin:prop%all_end)
        integer(int64)                         :: e_i
        integer                                :: i
        e_i = 0_int64
        !$acc parallel loop reduction(+:e_i)
        do i = prop%begin, prop%end
            e_i = e_i - ising_arr(i)*&
                ( ising_arr(i+1)&
                + ising_arr(i+prop%nx) )
        end do
        e = e_i / real(prop%particles, rkind)
        return
    end function ising2d_GPU_calc_energy
    !! ising2d_GPU_calc_magne: Calculate magnetism per site.
    pure real(rkind) function ising2d_GPU_calc_magne(prop, ising_arr) result(m)
        type(IG2_property)        , intent(in) :: prop
        integer           , device, intent(in) :: ising_arr(prop%all_begin:prop%all_end)
        integer(int64)                         :: m_i
        integer                                :: i
        m_i = 0_int64
        !$acc parallel loop reduction(+:m_i)
        do i = prop%begin, prop%end
            m_i = m_i + ising_arr(i)
        end do
        m = m_i / real(prop%particles, rkind)
        return
    end function ising2d_GPU_calc_magne
end module ising2d_GPU_m
