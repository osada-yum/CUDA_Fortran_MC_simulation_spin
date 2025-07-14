set -x -u -e

output_dir="data/xy2d"
outputfile="${output_dir}/xy2d_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_${execname}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "${output_dir}"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i \
    -e "s/nx = [0-9]*_int64/nx = ${nx}_int64/" \
    -e "s/mcs = [0-9]*_int32/mcs = ${mcs}_int32/" \
    -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
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
echo "xy2d_MET_GPU(${execname}),  ${nx}x${ny} == $((nx*ny)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a gpu_xy2d.log

rm -rf "${root_dir:?}"
