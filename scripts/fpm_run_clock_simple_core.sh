set -x -u -e

root_dir="./bin/root_$$"

## どのプログラムを使うかを選択する.
case ${kind} in
    "simple")
        srcfile="./src/clock/clock_simple_gpu_m.f90"
        progfile="./app/clock_simple_gpu_relaxation.f90"
        execname="clock_simple_gpu_relaxation"
        execfile="${root_dir}/bin/${execname}"
        ;;
    "table")
        srcfile="./src/clock/clock_table_gpu_m.f90"
        progfile="./app/clock_table_gpu_relaxation.f90"
        execname="clock_table_gpu_relaxation"
        execfile="${root_dir}/bin/${execname}"
        ;;
    "tableall")
        srcfile="./src/clock/clock_tableall_gpu_m.f90"
        progfile="./app/clock_tableall_gpu_relaxation.f90"
        execname="clock_tableall_gpu_relaxation"
        execfile="${root_dir}/bin/${execname}"
        ;;
    "dual_lattice_yhalf_tableall")
        srcfile="./src/clock/clock_dual_lattice_yhalf_tableall_m.f90"
        progfile="./app/clock_dual_lattice_yhalf_tableall_gpu_relaxation.f90"
        execname="clock_dual_lattice_yhalf_tableall_gpu_relaxation"
        execfile="${root_dir}/bin/${execname}"
        ;;
    *)
        if [ -z "${kind}" ]; then
            echo "Unknown kind: '${kind}'" 1>&2
        fi
        cat << EOF 1>&2
Selectable values:
simple
table
tableall
dual_lattice_yhalf_tableall
EOF
        exit 1
        ;;
esac

output_dir="data/${mstate}clock_periodic"
outputfile="${output_dir}/${mstate}clock_periodic_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_${execname}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "${output_dir}"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

## lock ファイルで同時にビルドを防ぐ.
lock_file="fpm_run.lock"
cnts=0
while [ -f "${lock_file}" ]; do
    sleep 30
    cnts=$((cnts+1))
    if echo "${cnts} > 100" | bc ; then
        break
    fi
done
touch "${lock_file}"
trap "{
     rm -f ${lock_file};
}" SIGINT SIGUSR1 SIGQUIT SIGKILL SIGSEGV


## プログラムを変更したり, ビルドしたりするのでlockが必要
## プログラムの定数の変更
sed -i -e "s/nx = [0-9]*_int64/nx = ${nx}_int64/" \
    -e "s/mstate = [0-9]*/mstate = ${mstate}/" \
    -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
    "${srcfile}"

sed -i -e "s/mcs = [0-9]*_int32/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = [0-9]*_int32/tot_sample = ${sample}_int32/" \
    -e "s/n_skip = [0-9]*_int64/n_skip = ${n_skip}_int64/" \
    -e "s/iseed = [0-9]*_int32/iseed = ${iseed}_int32/" \
    "${progfile}"

## ビルド
fpm clean --skip
fpm install "${execname}" --prefix="${root_dir}" --verbose --compiler='nvfortran' --flag="${FCFLAGS}"

## ビルドが終わったので, lockを解除する
kill -10 $$ # SIGUSR1
trap SIGINT
trap SIGUSR1
trap SIGQUIT
trap SIGKILL
trap SIGSEGV

## プログラム実行
start=$(date +%s)
${execfile} > "${tmpfile}"
end=$(date +%s)
cp -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 400 "${outputfile}"

## 経過時間の計算
hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
## 経過時間の計算
## logファイルへの書き込み.
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "${mstate}clock_MET_GPU(${execname}),  ${nx}x${ny} == $((nx*ny)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a "gpu_${mstate}clock.log"

## バイナリファイルの削除(最悪しなくても良い)
rm -rf "${root_dir:?}"
