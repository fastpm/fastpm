. testfunctions.sh

mpirun -n 4 $FASTPM cola.lua || fail

test -f cola/powerspec_1.0000.txt || fail
test -d cola/fastpm_1.0000 || fail
