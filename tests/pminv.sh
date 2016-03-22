. testfunctions.sh

mpirun -n 4 $FASTPM pminv.lua || fail
test -f pminv/powerspec_1.0000.txt || fail
test -d pminv/fastpm_1.0000 || fail
