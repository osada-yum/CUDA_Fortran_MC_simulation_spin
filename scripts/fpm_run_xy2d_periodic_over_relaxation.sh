set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=10000
ny=${nx}
mcs=10000
sample=500
kbt=0.895
iseed=42
n_skip=0

mcs_over_relax=${mcs}
n_over_relax=1

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_gpu_m.f90"
progfile="./app/xy2d_periodic_gpu_over_relaxation.f90"
execname="xy2d_periodic_gpu_over_relaxation"
execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_over_relaxation_core.sh"
