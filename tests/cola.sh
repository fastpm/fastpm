. testfunctions.sh

mpirun -n 4 $FASTPM cola.lua || fail

test -f cola/powerspec_1.0000.txt || fail
test -f cola/snp_1.0000.bin.00 || fail
