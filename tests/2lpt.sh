. testfunctions.sh

mpirun -n 4 $FASTPM 2lpt.lua || fail

test -f 2lpt/powerspec_1.0000.txt || fail
test -d 2lpt/fastpm_1.0000 || fail
