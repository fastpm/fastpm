#! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"
log=`mktemp`
check=`mktemp`

assert_success "mpirun -n 4 $FASTPM lightcone.lua > $log"

cat > $check << EOF
CHECK: Input power spectrum sigma8 0.815897

CHECK: Writing a catalog to lightcone/usmesh [1]
CHECK: Writing 0 objects.
CHECK: Writing a catalog to lightcone/usmesh [LL-0.200]
CHECK: Writing 0 objects.
CHECK: Total number of particles wrote into lightcone: 0

CHECK: Appending usmesh catalog to lightcone/usmesh
CHECK: Appending a catalog to lightcone/usmesh [1]
CHECK: Writing 423640 objects.
CHECK: Appending a catalog to lightcone/usmesh [LL-0.200]
CHECK: Writing 0 objects.

CHECK: Appending usmesh catalog to lightcone/usmesh
CHECK: Appending a catalog to lightcone/usmesh [1]
CHECK: Writing 607701 objects.
CHECK: Appending a catalog to lightcone/usmesh [LL-0.200]
CHECK: Writing 0 objects.
CHECK: Total number of particles wrote into lightcone: 1031341

CHECK: Appending usmesh catalog to lightcone/usmesh
CHECK: Appending a catalog to lightcone/usmesh [1]
CHECK: Writing 1523313 objects.
CHECK: Appending a catalog to lightcone/usmesh [LL-0.200]
CHECK: Writing 453 objects.

CHECK: Writing a catalog to lightcone/fastpm_1.0000 [1]
CHECK: Writing 2097152 objects.
CHECK: Writing a catalog to lightcone/fof_1.0000 [LL-0.200]
CHECK: Writing 22062 objects.

CHECK: Appending usmesh catalog to lightcone/usmesh
CHECK: Appending a catalog to lightcone/usmesh [1]
CHECK: Writing 1507785 objects.
CHECK: Appending a catalog to lightcone/usmesh [LL-0.200]
CHECK: Writing 4983 objects.
EOF

echo "---- Validating the log output $log with $check-------"
assert_success "cat $log | filecheck $check"
