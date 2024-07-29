set -x -u -e

FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=451
ny=${nx}
nz=$((nx - 1))
mcs=10000
sample=4700
kbt=4.51152174982078
iseed=42
n_skip=0

progfile="./app/ising3d_gpu_relaxation.f90"
execfile="./bin/ising3d_gpu_relaxation"
outputfile="data/3D-Ising/3D-Ising_GPU_x${nx}_y${ny}_z${nz}_mcs${mcs}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "data/3D-Ising"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/nz = .*/nz = ${nz}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    -e "s/iseed = .*/iseed = ${iseed}_int32/" \
    -e "s/n_skip = .*/n_skip = ${n_skip}_int64/" \
    "${progfile}"

fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" build ising3d_gpu_relaxation
start=$(date +%s)
time fpm install ising3d_gpu_relaxation --prefix="." --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
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
echo "3D-Ising_MET_GPU,  ${nx}x${ny}x${nz} == $((nx*ny*nz)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_ising3d.log
