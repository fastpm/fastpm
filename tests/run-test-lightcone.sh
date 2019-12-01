#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM lightcone.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Total number of particles wrote into lightcone: 3083305'
assert_file_contains $log 'Total number of particles wrote into lightcone: 7666753'
assert_file_contains $log 'Total number of particles wrote into lightcone: 9900645'
assert_file_contains $log 'Total number of particles wrote into lightcone: 11690327'
assert_file_contains $log 'Total number of particles wrote into lightcone: 12645089'
assert_file_contains $log 'Total number of particles wrote into lightcone: 13354288'
assert_file_contains $log 'Total number of particles wrote into lightcone: 13735853'
assert_file_contains $log 'Total number of particles wrote into lightcone: 13987789'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14119749'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14202578'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14236249'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14248346'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14251970'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14252293'

assert_file_contains $log 'Writing 7666753 objects.'
assert_file_contains $log 'Writing 4978336 objects.'
assert_file_contains $log 'Writing 453 objects.'
assert_file_contains $log 'Writing 2097152 objects.'
assert_file_contains $log 'Writing 22062 objects.'
assert_file_contains $log 'Writing 1607204 objects.'
assert_file_contains $log 'Writing 4983 objects.'
assert_file_contains $log 'sigma8.*0.815897'

report_test_status
