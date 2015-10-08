/* libfastpm: */
#include <mpi.h>

#include "libfastpm.h"
#include "pmtimer.h"
#include "msg.h"

int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm) 
{
    msg_init(comm);
    msg_set_loglevel(verbose);
    timer_init();

    pm_store_init(p);

    pm_store_alloc_evenly(p, pow(nc, 3), 2.0, comm);

    return 0;
}
