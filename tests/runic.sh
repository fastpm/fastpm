. testfunctions.sh
mpirun -n 4 $FASTPM ic.lua || fail
mpirun -n 4 $FASTPM runic.lua || fail
mpirun -n 4 $FASTPM runicnoise.lua || fail
mpirun -n 4 $FASTPM runicnoisek.lua || fail
test -f runic/powerspec00100_1.0000.txt || fail
test -f runic/snp00100_1.0000.bin.00 || fail
test -f runicnoise/powerspec00100_1.0000.txt || fail
test -f runicnoise/snp00100_1.0000.bin.00 || fail
test -f runicnoisek/powerspec00100_1.0000.txt || fail
test -f runicnoisek/snp00100_1.0000.bin.00 || fail
