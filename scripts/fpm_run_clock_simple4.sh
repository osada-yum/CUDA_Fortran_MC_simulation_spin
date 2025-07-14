set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=2000
ny=${nx}
mcs=10000
sample=1000
kbt=0.450
iseed=42
n_skip=0
mstate=6

script_dir="scripts"
kind="dual_lattice_yhalf_tableall"
source "${script_dir}/fpm_run_clock_simple_core.sh"
