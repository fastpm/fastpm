. testfunctions.sh

mpirun -n 4 $FASTPM pm.lua || fail
test -f pm/powerspec00100_1.0000.txt || fail
test -f pm/snp00100_1.0000.bin.00 || fail
