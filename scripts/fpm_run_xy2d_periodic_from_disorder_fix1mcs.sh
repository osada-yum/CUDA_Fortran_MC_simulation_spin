set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=1500
ny=${nx}
mcs=100000
sample=2000 #12500
kbt=0.890
iseed=42
n_skip=0

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_gpu_m.f90"
progfile="./app/xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs.f90"
execname="xy2d_periodic_gpu_relaxation_from_disorder_fix1mcs"

execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_from_disorder_core.sh"
