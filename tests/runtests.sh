#! /bin/bash

source testfunctions.sh

set -x
mpirun -n 4 $FASTPM standard.lua ic || fail

mpirun -n 4 $FASTPM standard.lua fastpm || fail
mpirun -n 4 $FASTPM standard.lua fastpm lanczos3 || fail
exit
mpirun -n 4 $FASTPM standard.lua fastpm gaussian || fail
mpirun -n 4 $FASTPM standard.lua fastpm aggressive || fail
mpirun -n 4 $FASTPM standard.lua fastpm gadget || fail
mpirun -n 4 $FASTPM standard.lua fastpm eastwood || fail
mpirun -n 4 $FASTPM standard.lua fastpm eastwood gaussian|| fail
mpirun -n 4 $FASTPM standard.lua fastpm twothird || fail
mpirun -n 4 $FASTPM standard.lua fastpm gaussian36 || fail

mpirun -n 4 $FASTPM standard.lua za fnl|| fail

mpirun -n 4 $FASTPM standard.lua fastpm lineark || fail
mpirun -n 4 $FASTPM standard.lua fastpm whitenoisek || fail

mpirun -n 4 $FASTPM standard.lua za || fail
mpirun -n 4 $FASTPM standard.lua 2lpt || fail
mpirun -n 4 $FASTPM standard.lua fastpm || fail
mpirun -n 4 $FASTPM standard.lua pm || fail
mpirun -n 4 $FASTPM standard.lua cola || fail

mpirun -n 4 $FASTPM standard.lua fastpm inverted || fail
mpirun -n 4 $FASTPM standard.lua fastpm remove_variance || fail

