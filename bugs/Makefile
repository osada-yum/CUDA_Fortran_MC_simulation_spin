# Bug occur
# nvfortran version 25.1 (Ubuntu 20.04.6 LTS, from apt "nvhpc-25-1")
FC = /opt/nvidia/hpc_sdk/Linux_x86_64/25.1/compilers/bin/nvfortran
FCFLAGS = -g -Mchkptr -traceback -Wall -Wextra -acc -cuda -cudalib=curand -Minfo=accel

all: bug_xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs.out

bug_xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs.out: bug_xy2d_periodic_gpu_m.o bug_xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs.f90
	$(FC) $(FCFLAGS) $^ -o $@

bug_xy2d_periodic_gpu_m.o: bug_xy2d_periodic_gpu_m.f90
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.mod

distclean:
	rm -rf *.o *.mod *.out
