set -x -u -e

root_dir="./bin/root_$$"

# srcfile="./src/clock_simple_gpu_m.f90"
# progfile="./app/clock_simple_gpu_relaxation.f90"
# execname="clock_simple_gpu_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/clock_table_gpu_m.f90"
# progfile="./app/clock_table_gpu_relaxation.f90"
# execname="clock_table_gpu_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/clock_tableall_gpu_m.f90"
# progfile="./app/clock_tableall_gpu_relaxation.f90"
# execname="clock_tableall_gpu_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/clock_dual_lattice_tableall_m.f90"
# progfile="./app/clock_dual_lattice_tableall_gpu_relaxation.f90"
# execname="clock_dual_lattice_tableall_gpu_relaxation"
# execfile="${root_dir}/bin/${execname}"

srcfile="./src/clock_dual_lattice_yhalf_tableall_m.f90"
progfile="./app/clock_dual_lattice_yhalf_tableall_gpu_relaxation.f90"
execname="clock_dual_lattice_yhalf_tableall_gpu_relaxation"
execfile="${root_dir}/bin/${execname}"

output_dir="data/6clock_periodic"
outputfile="${output_dir}/6clock_periodic_test_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_${execname}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "${output_dir}"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = [0-9]*_int64/nx = ${nx}_int64/" \
    -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
    "${srcfile}"

sed -i -e "s/mcs = [0-9]*_int32/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = [0-9]*_int32/tot_sample = ${sample}_int32/" \
    -e "s/n_skip = [0-9]*_int64/n_skip = ${n_skip}_int64/" \
    -e "s/iseed = [0-9]*_int32/iseed = ${iseed}_int32/" \
    "${progfile}"

fpm install "${execname}" --prefix="${root_dir}" --verbose --compiler='nvfortran' --flag="${FCFLAGS}"
start=$(date +%s)
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
echo "6clock_MET_GPU(${execname}, test),  ${nx}x${ny} == $((nx*ny)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_6clock.log

rm -rf "${root_dir:?}"
