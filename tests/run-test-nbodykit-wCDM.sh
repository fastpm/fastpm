#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM nbodykit-wCDM.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Writing 1918 objects.'
assert_file_contains $log 'Writing 1506 objects.'
assert_file_contains $log 'RSD factor.*1.162687e-02'
assert_file_contains $log 'sigma8.*0.815897'

report_test_status
