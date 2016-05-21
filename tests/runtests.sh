#! /bin/bash

source testfunctions.sh

set -x
mpirun -n 4 $FASTPM standard.lua ic || fail
mpirun -n 4 $FASTPM standard.lua pm lineark || fail
mpirun -n 4 $FASTPM standard.lua pm whitenoisek || fail

mpirun -n 4 $FASTPM standard.lua za || fail
mpirun -n 4 $FASTPM standard.lua 2lpt || fail
mpirun -n 4 $FASTPM standard.lua pm || fail
mpirun -n 4 $FASTPM standard.lua zola || fail
mpirun -n 4 $FASTPM standard.lua cola || fail

mpirun -n 4 $FASTPM standard.lua zola inverted || fail
mpirun -n 4 $FASTPM standard.lua zola remove_variance || fail

