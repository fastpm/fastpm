#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM restart.lua > $log"
assert_success "mpirun -n 4 $FASTPM restart.lua -r restart/fastpm_0.5000> $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Velocity dispersion (a = 0.6124): std = 1.63807 1.75754 1.94999'
assert_file_contains $log 'Velocity dispersion (a = 0.8660): std = 2.44703 2.62561 2.90857'

echo "Testing conversion from Gadget"
assert_success "python ../python/convert-to-gadget-1.py restart/fastpm_0.5000 restart/gadget_0.5000"
assert_success "python ../python/convert-from-gadget-1.py 'restart/gadget_0.5000.*' restart/fastpm_0.5000_from_gadget"

assert_success "mpirun -n 4 $FASTPM restart.lua -r restart/fastpm_0.5000_from_gadget> $log"
echo "---- Validating the log output -------"
assert_file_contains $log 'Velocity dispersion (a = 0.6124): std = 1.63807 1.75754 1.94999'
assert_file_contains $log 'Velocity dispersion (a = 0.8660): std = 2.44703 2.62561 2.90857'

report_test_status
