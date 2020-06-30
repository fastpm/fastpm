#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`

assert_success "mpirun -n 4 $FASTPM lightcone-ODE.lua > $log"

echo "---- Validating the log output -------"
assert_file_contains $log 'Total number of particles wrote into lightcone: 5540104'
assert_file_contains $log 'Total number of particles wrote into lightcone: 10500104'
assert_file_contains $log 'Total number of particles wrote into lightcone: 12733996'
assert_file_contains $log 'Total number of particles wrote into lightcone: 14523678'
assert_file_contains $log 'Total number of particles wrote into lightcone: 15478440'
assert_file_contains $log 'Total number of particles wrote into lightcone: 16187639'
assert_file_contains $log 'Total number of particles wrote into lightcone: 16569204'
assert_file_contains $log 'Total number of particles wrote into lightcone: 16821140'
assert_file_contains $log 'Total number of particles wrote into lightcone: 16953100'
assert_file_contains $log 'Total number of particles wrote into lightcone: 17035929'
assert_file_contains $log 'Total number of particles wrote into lightcone: 17069600'
assert_file_contains $log 'Total number of particles wrote into lightcone: 17081697'
assert_file_contains $log 'Total number of particles wrote into lightcone: 17085321'
assert_file_contains $log 'Total number of particles wrote into lightcone: 17085644'

assert_file_contains $log 'Writing 5540104 objects.'
assert_file_contains $log 'Writing 4960000 objects.'
assert_file_contains $log 'Writing 4978336 objects.'
assert_file_contains $log 'Writing 453 objects.'
assert_file_contains $log 'Writing 2097152 objects.'
assert_file_contains $log 'Writing 22063 objects.'
assert_file_contains $log 'Writing 1607204 objects.'
assert_file_contains $log 'Writing 4983 objects.'

assert_file_contains $log 'sigma8.*0.815897'

report_test_status
