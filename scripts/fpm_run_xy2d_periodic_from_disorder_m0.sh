set -x -u -e

# m0は上手くいかない... 磁化 m(t) / m(0) がガタつくし, 非単調な振る舞い.

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=2000
ny=${nx}
mcs=10000
sample=12500
kbt=0.890
iseed=42
n_skip=0

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_gpu_m.f90"
progfile="./app/xy2d_periodic_gpu_relaxation_from_disorder_m0.f90"
execname="xy2d_periodic_gpu_relaxation_from_disorder_m0"

execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
source "${script_dir}/fpm_run_xy2d_periodic_from_disorder_core.sh"
