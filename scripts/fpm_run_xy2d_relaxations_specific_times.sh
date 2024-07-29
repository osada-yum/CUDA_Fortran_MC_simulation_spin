#!/bin/bash

set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=1001
ny=${nx}
mcs=100000
bin=1000
sample=1
kbt=0.90
iseed=42
n_skip=0

root_dir="/tmp"
progfile="./app/xy2d_gpu_relaxations_specific_times.f90"
execfile="${root_dir}/bin/xy2d_gpu_relaxations_specific_times"
outputdir="data/2D-XY_relaxations"
outputfile="${outputdir}/2D-XY_GPU_x${nx}_y${ny}_specific_mcs${mcs}_bin${bin}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p ${outputdir}
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/bin = .*/bin = ${bin}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    -e "s/iseed = .*/iseed = ${iseed}_int32/" \
    -e "s/n_skip = .*/n_skip = ${n_skip}_int64/" \
    "${progfile}"

fpm clean --skip
time fpm install xy2d_gpu_relaxations --prefix="${root_dir}" --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
start=$(date +%s)
${execfile} > "${tmpfile}"
end=$(date +%s)
mv -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 400 "${outputfile}"

hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "2D-XY_MET_GPU_relaxations_specific_times,  ${nx}x${ny} == $((nx*ny)), bin${bin}_sample${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_xy2d_relaxations.log
