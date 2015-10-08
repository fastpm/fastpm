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

int fastpm_particle_to_mesh(PM * pm, PMStore * p) {
    /* After this function, pm->canvas contains the real space density field */
    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
    };

    pm_append_ghosts(&pgd);

    pm_paint(pm, p, p->np + pgd.nghosts);

    pm_destroy_ghosts(&pgd);
}
