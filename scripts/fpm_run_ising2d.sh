set -x -u -e

FCFLAGS="-O2 -acc -cuda -gpu=cc89 -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
progfile="./app/ising2d_gpu_relaxation.f90"
nx=10001
ny=$((nx - 1))
mcs=1000
sample=1
kbt=2.269185314213022
iseed=42

progfile="./app/ising2d_gpu_relaxation.f90"
execfile="./bin/ising2d_gpu_relaxation"
outputfile="data/2D-Ising/2D-Ising_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "data/2D-Ising"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/iseed = .*/iseed = ${iseed}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    "${progfile}"

start=$(date +%s)
time fpm install ising2d_gpu_relaxation --prefix="." --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
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
echo "2D-Ising_MET_GPU,  ${nx}x${ny} == $((nx*ny)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_ising2d.log
