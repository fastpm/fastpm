. testfunctions.sh

mpirun -n 4 $FASTPM za.lua || fail

test -f za/powerspec_1.0000.txt || fail
test -f za/snp_1.0000.bin.00 || fail
