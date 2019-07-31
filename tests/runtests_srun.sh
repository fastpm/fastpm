#! /bin/bash

source testfunctions.sh

set -x
srun -n 4 $FASTPM standard.lua za || die
srun -n 4 $FASTPM standard.lua 2lpt || die
srun -n 4 $FASTPM standard.lua pm || die
srun -n 4 $FASTPM standard.lua zola || die
srun -n 4 $FASTPM standard.lua cola || die

srun -n 4 $FASTPM standard.lua zola inverted || die
srun -n 4 $FASTPM standard.lua zola remove_variance || die

