set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=10001
ny=$((nx - 1))
mcs=10000
sample=1
kbt=0.895
iseed=42
n_skip=0

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_core.sh"
