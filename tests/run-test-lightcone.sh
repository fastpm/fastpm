#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM lightcone.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Total number of particles in light cone slice: 3083305'
assert_file_contains $log 'Total number of particles in light cone slice: 2233892'
assert_file_contains $log 'Total number of particles in light cone slice: 4023574'
assert_file_contains $log 'Total number of particles in light cone slice: 709199'
assert_file_contains $log 'Total number of particles in light cone slice: 1090764'
assert_file_contains $log 'Total number of particles in light cone slice: 1342700'
assert_file_contains $log 'Total number of particles in light cone slice: 1474660'
assert_file_contains $log 'Total number of particles in light cone slice: 1557489'
assert_file_contains $log 'Total number of particles in light cone slice: 1591160'
assert_file_contains $log 'Total number of particles in light cone slice: 1603257'
assert_file_contains $log 'Total number of particles in light cone slice: 1607204'

assert_file_contains $log 'Writing 7666753 objects.'
assert_file_contains $log 'Writing 4978336 objects.'
assert_file_contains $log 'Writing 453 objects.'
assert_file_contains $log 'Writing 2097152 objects.'
assert_file_contains $log 'Writing 22062 objects.'
assert_file_contains $log 'Writing 1607204 objects.'
assert_file_contains $log 'Writing 4983 objects.'
assert_file_contains $log 'sigma8.*0.815897'

report_test_status
