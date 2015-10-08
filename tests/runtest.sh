export OMP_NUM_THREADS=1
mkdir -p example-output
mpirun -n 4 ../fastpm ./example.lua
