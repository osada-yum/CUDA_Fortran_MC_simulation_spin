set -x -u -e

FCFLAGS="-O3 -use_fast_math -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=1001
ny=$((nx - 1))
mcs=10000
sample=5000
kbt=0.95

progfile="./app/xy2d_gpu_relaxation.f90"
execfile="./bin/xy2d_gpu_relaxation"
outputfile="data/XY/2D-XY_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_$(date +%Y%m%d_%H%M%S).dat"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    "${progfile}"

fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" build xy2d_gpu_relaxation
start=$(date +%s)
time fpm install xy2d_gpu_relaxation --prefix="." --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
${execfile} > "${tmpfile}"
end=$(date +%s)
cp -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 400 "${outputfile}"

hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "2D-XY_MET_GPU,  ${nx}x${ny} == $((nx*ny)),»  ${sample},   ${mcs},»   ${kbt}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_xy2d.log
