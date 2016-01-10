. testfunctions.sh

mpirun -n 4 $FASTPM 2lpt|| fail

test -f 2lpt/powerspec00100_1.0000.txt || fail
test -f 2lpt/snp00100_1.0000.bin.00 || fail
