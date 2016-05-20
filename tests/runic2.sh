. testfunctions.sh
srun -n 4 $FASTPM ic.lua || fail
srun -n 4 $FASTPM runic.lua || fail
test -f runic/powerspec_1.0000.txt || fail
test -d runic/fastpm_1.0000 || fail
srun -n 4 $FASTPM runicnoise.lua || fail
test -f runicnoise/powerspec_1.0000.txt || fail
test -d runicnoise/fastpm_1.0000 || fail
srun -n 4 $FASTPM runicnoisek.lua || fail
test -f runicnoisek/powerspec_1.0000.txt || fail
test -d runicnoisek/fastpm_1.0000 || fail
