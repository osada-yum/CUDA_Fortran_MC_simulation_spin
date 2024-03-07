set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=501
ny=$((nx - 1))
mcs=100000
sample=150
kbt=0.80
state=6
n_multi=2

progfile="./app/clock_gpu_multi_relaxation.f90"
execfile="./bin/clock_gpu_multi_relaxation"
outputfile="data/clock/${state}state-clock_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}x${n_multi}_kbt${kbt}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "data/clock"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    -e "s/state = .*/state = ${state}_int32/" \
    -e "s/n_multi = .*/n_multi = ${n_multi}_int32/" \
    "${progfile}"

fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" build clock_gpu_relaxation
start=$(date +%s)
time fpm install clock_gpu_relaxation --prefix="." --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
${execfile} > "${tmpfile}"
end=$(date +%s)
cp -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 600 "${outputfile}"

hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "${state}state-clock_MET_GPU,  ${nx}x${ny} == $((nx*ny)),»  ${sample}x${n_multi} == $((sample*n_multi)),   ${mcs},»   ${kbt}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_clock.log
