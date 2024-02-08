FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
progfile="./app/xy2d_gpu_relaxation.f90"
nx=1001
ny=$((nx - 1))
mcs=1000
sample=1440000
kbt=1.0
outputfile="data/XY/2D-XY_GPU_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_$(date +%Y%m%d_%H%M%S).dat"
sed -i -e "s/nx = .*/nx = ${nx}_int64/" \
    -e "s/ny = .*/ny = ${ny}_int64/" \
    -e "s/mcs = .*/mcs = ${mcs}_int32/" \
    -e "s/tot_sample = .*/tot_sample = ${sample}_int32/" \
    -e "s/kbt = .*/kbt = ${kbt}_real64/" \
    "${progfile}"

fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" build xy2d_gpu_relaxation
fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" run xy2d_gpu_relaxation | tee ${outputfile}
