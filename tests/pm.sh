. testfunctions.sh

mpirun -n 4 $FASTPM pm.lua || fail
test -f pm/powerspec_1.0000.txt || fail
test -d pm/fastpm_1.0000 || fail
