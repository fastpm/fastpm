#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM rfof.lua > $log"

echo "---- Validating the log output -------"
# FIXME: Update these numbers after checking if the model is sane.
assert_file_contains $log 'Writing 7839 objects.'
assert_file_contains $log 'Writing 15165 objects.'
assert_file_contains $log 'RSD factor.*1.140331e-02'
assert_file_contains $log 'sigma8.*0.815897'

report_test_status
