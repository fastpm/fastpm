export OMP_NUM_THREADS=1
mpirun -n 4 ../src/fastpm ./example.lua
