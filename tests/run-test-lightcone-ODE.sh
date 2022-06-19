#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`
check=run-test-lightcone-ODE.check

assert_success "mpirun -n 4 $FASTPM lightcone-ODE.lua > $log"

echo "---- Validating the log output $log with $check-------"
assert_success "cat $log | filecheck $check"
