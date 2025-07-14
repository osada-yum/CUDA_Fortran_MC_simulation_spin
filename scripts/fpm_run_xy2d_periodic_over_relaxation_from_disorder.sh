#!/bin/bash

set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -Minfo=accel"
# FCFLAGS="${FCFLAGS} -g -Mbounds -gpu=debug"
nx=1500
ny=${nx}
mcs=1000
sample=40000
kbt=0.890
iseed=42
n_skips=(40 41 42 43 44 45)

n_over_relax=1

root_dir="./bin/root_$$"

srcfile="./src/xy2d_periodic_gpu_m.f90"
# progfile="./app/xy2d_periodic_gpu_over_relaxation.f90"
# execname="xy2d_periodic_gpu_over_relaxation"

progfile="./app/xy2d_periodic_gpu_over_relaxation_from_disorder.f90"
execname="xy2d_periodic_gpu_over_relaxation_from_disorder"

execfile="${root_dir}/bin/${execname}"

script_dir="scripts"
for n_skip in ${n_skips[@]}
do
    source "${script_dir}/fpm_run_xy2d_periodic_over_relaxation_from_disorder_core.sh"
done
