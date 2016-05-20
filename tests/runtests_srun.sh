#! /bin/bash

source testfunctions.sh

set -x
srun -n 4 $FASTPM standard.lua za || fail
srun -n 4 $FASTPM standard.lua 2lpt || fail
srun -n 4 $FASTPM standard.lua pm || fail
srun -n 4 $FASTPM standard.lua zola || fail
srun -n 4 $FASTPM standard.lua cola || fail

srun -n 4 $FASTPM standard.lua zola inverted || fail
srun -n 4 $FASTPM standard.lua zola remove_variance || fail

