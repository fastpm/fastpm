#assert_success "! /bin/bash

source testfunctions.sh

FASTPM="`dirname $0`/../src/fastpm -T 1"

assert_success "mpirun -n 4 $FASTPM standard.lua ic"
assert_success "mpirun -n 4 $FASTPM standard.lua fastpm lineark"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm whitenoisek"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm"

#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm lanczos2"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm lanczos3"

#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm gaussian"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm aggressive"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm gadget"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm eastwood"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm eastwood gaussian"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm twothird"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm gaussian36"

assert_success "mpirun -n 4 $FASTPM standard.lua za fnl"


assert_success "mpirun -n 4 $FASTPM standard.lua za"
#assert_success "mpirun -n 4 $FASTPM standard.lua 2lpt"
assert_success "mpirun -n 4 $FASTPM standard.lua fastpm"
#assert_success "mpirun -n 4 $FASTPM standard.lua pm"
#assert_success "mpirun -n 4 $FASTPM standard.lua cola"

#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm inverted"
#assert_success "mpirun -n 4 $FASTPM standard.lua fastpm remove_variance"

