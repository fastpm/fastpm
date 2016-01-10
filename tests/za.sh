. testfunctions.sh

mpirun -n 4 $FASTPM za|| fail

test -f za/powerspec00100_1.0000.txt || fail
test -f za/snp00100_1.0000.bin.00 || fail
