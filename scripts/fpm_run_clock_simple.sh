set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=1000
ny=${nx}
mcs=100000
sample=250
kbt=0.91
iseed=42
n_skip=0

script_dir="scripts"
source "${script_dir}/fpm_run_clock_simple_core.sh"
