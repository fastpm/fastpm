#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

log=lclog
assert_success "mpirun -n 4 $FASTPM lightcone.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Writing 1607204 objects.'
assert_file_contains $log 'Writing 4983 objects.'
assert_file_contains $log 'sigma8.*0.815897'

report_test_status
