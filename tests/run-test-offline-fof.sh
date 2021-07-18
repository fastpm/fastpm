#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
FOF="`dirname $0`/../src/fastpm-fof -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM restart.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Velocity dispersion (a = 0.6124): std = 1.63807 1.75754 1.94999'
assert_file_contains $log 'Velocity dispersion (a = 0.8660): std = 2.44703 2.62561 2.90857'
assert_file_contains $log 'Writing 4017 objects'

log=`mktemp`
assert_success "mpirun -n 4 $FOF restart/fastpm_1.0000 0.2 > $log"
echo "---- Validating the log output -------"
assert_file_contains $log 'Writing 4017 objects'

report_test_status
