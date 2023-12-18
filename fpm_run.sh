FCFLAGS="-O2 -acc -cuda -cudalib=curand"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
fpm --verbose --compiler='nvfortran' --flag="${FCFLAGS}" run
