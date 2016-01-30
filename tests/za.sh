. testfunctions.sh

mpirun -n 4 $FASTPM za.lua || fail

test -f za/powerspec_1.0000.txt || fail
test -d za/fastpm_1.0000 || fail
