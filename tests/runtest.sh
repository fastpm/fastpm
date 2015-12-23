export OMP_NUM_THREADS=1
fail () {
    exit 1
}
mpirun -n 4 ../src/fastpm cola.lua || fail
mpirun -n 4 ../src/fastpm pm.lua || fail
