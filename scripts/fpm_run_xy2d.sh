set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=10001
ny=$((nx - 1))
mcs=10000
sample=500
kbt=0.890
iseed=42
n_skip=0

root_dir="./bin/root_$$"

srcfile="./src/xy2d_gpu_m.f90"
progfile="./app/xy2d_gpu_relaxation.f90"
execname="xy2d_gpu_relaxation"
execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_core.sh"
