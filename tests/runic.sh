. testfunctions.sh
mpirun -n 4 $FASTPM ic.lua || fail
mpirun -n 4 $FASTPM runic.lua || fail
test -f runic/powerspec00100_1.0000.txt || fail
test -f runic/snp00100_1.0000.bin.00 || fail
