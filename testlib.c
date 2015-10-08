#include <stdio.h>
#include <string.h>
#include "libfastpm.h"

void run2lpt(PM * pm, PMStore * pdata, double a, double omega_m, 
        real_t * deltak_0, real_t * deltak_1, MPI_Comm comm) {

    pm_start(pm);

    memcpy(pm->canvas, deltak_0, sizeof(pm->canvas[0]) * pm->allocsize);

    pm_2lpt_main(pm, pdata, comm);

    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */
    double shift[3] = {0, 0, 0};

    pm_2lpt_set_initial(a, omega_m, pdata, shift);

    fastpm_particle_to_mesh(pm, pdata);
    pm_r2c(pm);

    memcpy(deltak_1, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);
}

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

    real_t * deltak_0 = malloc(sizeof(deltak_0[0]) * pm->allocsize);
    real_t * deltak_1 = malloc(sizeof(deltak_1[0]) * pm->allocsize);

    int mcmcstep;

    for(mcmcstep = 0; mcmcstep < 1; mcmcstep++) {
        /* TODO: fill in the modes to canvas */
        memset(deltak_0, 0, sizeof(deltak_0[0]) * pm->allocsize);
        run2lpt(pm, pdata, 1.0, omega_m, deltak_0, deltak_1, comm);
    }
    free(deltak_0);
    free(deltak_1);

    pm_destroy(pm);

    MPI_Finalize();
    return 0;
}

