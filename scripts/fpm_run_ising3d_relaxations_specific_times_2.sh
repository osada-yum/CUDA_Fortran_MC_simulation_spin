#!/bin/bash

set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=451
ny=${nx}
nz=$((nx - 1))
mcs=10000
bin=1000
sample=1
kbt=4.51152174982078
iseed=42
n_skip=0

dest_install="/tmp"
progfile="./app/ising3d_gpu_relaxations_specific_times.f90"
execfile="${dest_install}/bin/ising3d_gpu_relaxations_specific_times"
outputdir="data/3D-Ising_relaxations"
outputfile="${outputdir}/3D-Ising_GPU_x${nx}_y${ny}_z${nz}_specific_mcs${mcs}_bin${bin}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p ${outputdir}
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/nz = .*/nz = ${nz}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/bin = .*/bin = ${bin}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    -e "s/iseed = .*/iseed = ${iseed}_int32/" \
    -e "s/n_skip = .*/n_skip = ${n_skip}_int64/" \
    "${progfile}"

fpm clean --skip
time fpm install ising3d_gpu_relaxations --prefix="${dest_install}" --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
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
echo "3D-Ising_MET_GPU_relaxations_specific_times,  ${nx}x${ny}x${nz} == $((nx*ny*nz)), bin${bin}_sample${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_ising3d_relaxations.log
