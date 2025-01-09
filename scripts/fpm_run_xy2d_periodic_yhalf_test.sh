set -x -u -e

# FCFLAGS="-Ofast -use_fast_math -acc -cuda -cudalib=curand"
FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -C -gpu=debug"
nx=10000
ny=${nx}
mcs=10000
sample=1
kbt=0.895
iseed=42
n_skip=0

root_dir="./bin/root_$$"

# srcfile="./src/xy2d_periodic_yhalf_gpu_m.f90"
# progfile="./app/xy2d_periodic_yhalf_gpu_relaxation.f90"
# execname="xy2d_periodic_yhalf_gpu_relaxation"

srcfile="./src/xy2d_periodic_yhalf_index_flip_gpu_m.f90"
progfile="./app/xy2d_periodic_yhalf_index_flip_gpu_relaxation.f90"
execname="xy2d_periodic_yhalf_index_flip_gpu_relaxation"

execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_core.sh"
