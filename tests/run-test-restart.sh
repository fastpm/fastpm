#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM restart.lua > $log"
assert_success "mpirun -n 4 $FASTPM restart.lua -r restart/fastpm_0.5000> $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Velocity dispersion (a = 0.6124): std = 1.63806 1.75753 1.94998'
assert_file_contains $log 'Velocity dispersion (a = 0.8660): std = 2.44702 2.6256 2.90856'

report_test_status
