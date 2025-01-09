set -x -u -e

echo "Wrong calculation?"

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=2000
ny=${nx}
mcs=100
sample=1
kbt=0.895
iseed=42
n_skip=0

NUM_THREADS=32

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_samples_gpu_m.f90"
progfile="./app/xy2d_periodic_samples_gpu_relaxation.f90"
execname="xy2d_periodic_samples_gpu_relaxation"
execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_samples_core.sh"
