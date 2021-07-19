#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
RFOF="`dirname $0`/../src/fastpm-rfof -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM restart.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Velocity dispersion (a = 0.6124): std = 1.63807 1.75754 1.94999'
assert_file_contains $log 'Velocity dispersion (a = 0.8660): std = 2.44703 2.62561 2.90857'
assert_file_contains $log 'Writing 14534 objects'

log=`mktemp`
assert_success "mpirun -n 4 $RFOF restart/fastpm_1.0000 1.0 > $log"
echo "---- Validating the log output -------"
assert_file_contains $log 'Writing 14533 objects'

report_test_status
