set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=1000
ny=${nx}
mcs=100
sample=500
kbt=0.890
iseed=42
n_skip=0
finite_magne=0.02d0

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_gpu_m.f90"
progfile="./app/xy2d_periodic_gpu_relaxation_from_disorder_finite_magne.f90"
execname="xy2d_periodic_gpu_relaxation_from_disorder_finite_magne"

execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_from_disorder_finite_magne_core.sh"
