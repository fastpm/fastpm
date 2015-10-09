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

    real_t * deltak_0 = malloc(sizeof(deltak_0[0]) * pm->allocsize);
    real_t * deltak_1 = malloc(sizeof(deltak_1[0]) * pm->allocsize);

    int mcmcstep;

    for(mcmcstep = 0; mcmcstep < 1; mcmcstep++) {
        /* Fill in some reasonable (EH PowerSpec) to the canvas */
        struct fastpm_powerspec_eh_params eh = {
            .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
            .hubble_param = 0.7,
            .omegam = 0.260,
            .omegab = 0.044,
        };
        fastpm_fill_deltak(pm, deltak_0, 100, fastpm_powerspec_eh, &eh);
        fastpm_evolve_2lpt(pm, pdata, 1.0, omega_m, deltak_0, deltak_1, comm);

        /* example looping over ks */
        ptrdiff_t ind;
        ptrdiff_t i[3] = {0};

        for(ind = 0; ind < pm->allocsize; ind += 2, pm_inc_o_index(pm, i)) {
            int d;
            double k[3];
            for(d = 0; d < 3; d ++) {
                k[d] = pm->MeshtoK[d][i[d] + pm->ORegion.start[d]];
            }
            deltak_1[ind + 0] = 0; /* real */
            deltak_1[ind + 1] = 0; /* imag */
        }
    }

    free(deltak_0);
    free(deltak_1);

    MPI_Finalize();
    return 0;
}

