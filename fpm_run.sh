FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
progfile="./app/ising2d_gpu_relaxation.f90"
nx=1001
ny=$((nx - 1))
mcs=1000
sample=1440000
kbt=2.26918531421
outputfile="data/2D-Ising_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_$(date +%Y%m%d_%H%M%S).dat"
sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    "${progfile}"

fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" build
fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" run | tee ${outputfile}
