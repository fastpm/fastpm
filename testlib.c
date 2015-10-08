#include <stdio.h>
#include <string.h>
#include "libfastpm.h"

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;
    int nc = 128;
    double boxsize = 256;
    double alloc_factor = 2.0;
    double omega_m = 0.292;

    PMStore * pdata = alloca(sizeof(PMStore));
    PM * pm = alloca(sizeof(PM));

    fastpm_init(pdata, nc, alloc_factor, comm);

    pm_2lpt_init(pm, pdata, nc, boxsize, comm);

    pm_start(pm);

    /* TODO: fill in the modes */

    memset(pm->canvas, 0, sizeof(pm->canvas[0]) * pm->allocsize);
    memset(pm->workspace, 0, sizeof(pm->workspace[0]) * pm->allocsize);

    pm_2lpt_main(pm, pdata, comm);

    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */
    double shift[3] = {0, 0, 0};

    pm_2lpt_set_initial(1.0, omega_m, pdata, shift);

    pm_destroy(pm);

    MPI_Finalize();
    return 0;
}

