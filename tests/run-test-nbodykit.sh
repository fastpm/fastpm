#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
mpirun -n 4 $FASTPM nbodykit.lua | tee log || fail
echo "---- Validating the output log file -------"
expect_contains log 'Writing 1894 objects.'
expect_contains log 'Writing 1668 objects.'
expect_contains log 'RSD factor.*1.140331e-02'
expect_contains log 'RSD factor.*1.140331e-02'
expect_contains log 'sigma8.*0.815897'

