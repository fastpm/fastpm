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

    fastpm_init_pm(pm, pdata, nc, boxsize, comm);

    real_t * rhok_0 = malloc(sizeof(rhok_0[0]) * pm->allocsize);
    real_t * rhok_1 = malloc(sizeof(rhok_1[0]) * pm->allocsize);
    real_t * rhok_1truth = malloc(sizeof(rhok_1truth[0]) * pm->allocsize);
    real_t * rhok_0truth = malloc(sizeof(rhok_0truth[0]) * pm->allocsize);
    real_t * Fk = malloc(sizeof(Fk[0]) * pm->allocsize);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_fill_deltak(pm, rhok_0truth, 100, (fastpm_pkfunc)fastpm_powerspec_eh, &eh);
    fastpm_evolve_2lpt(pm, pdata, 1.0, omega_m, rhok_0truth, rhok_1truth, comm);

    /* now a fake MCMC loop. */
    int mcmcstep;

    for(mcmcstep = 0; mcmcstep < 1; mcmcstep++) {
        /* example looping over ks */
        ptrdiff_t ind;
        ptrdiff_t i[3] = {0};
        double fac = sqrt(pow(2 * M_PI, 3) / pm->Volume);

        for(ind = 0; ind < pm->allocsize; ind += 2, pm_inc_o_index(pm, i)) {
            int d;
            double k[3];
            double kk = 0.0;
            for(d = 0; d < 3; d ++) {
                k[d] = pm->MeshtoK[d][i[d] + pm->ORegion.start[d]];
                kk += k[d] * k[d];
            }
            double k1 = sqrt(kk);
         
            double P = fastpm_powerspec_eh(sqrt(kk), &eh);
            rhok_0[ind + 0] = fac * sqrt(P) * 1.0; /* real */
            rhok_0[ind + 1] = fac * sqrt(P) * 1.0; /* imag */
        }

        /* Evolve to rhok_1 */
        fastpm_evolve_2lpt(pm, pdata, 1.0, omega_m, rhok_0, rhok_1, comm);

        /* The difference rhok_d */
        i[0] = 0;
        i[1] = 0;
        i[2] = 0;

        for(ind = 0; ind < pm->allocsize; ind += 2, pm_inc_o_index(pm, i)) {
            int d;
            double k[3];
            for(d = 0; d < 3; d ++) {
                k[d] = pm->MeshtoK[d][i[d] + pm->ORegion.start[d]];
            }
            rhok_1[ind + 0] -= rhok_1truth[ind + 0];
            rhok_1[ind + 1] -= rhok_1truth[ind + 1];
        }
        /* calculate the HMC force into Fk */
        fastpm_derivative_2lpt(pm, pdata, rhok_1, Fk, comm);
    }

    free(rhok_0);
    free(rhok_1);
    free(rhok_0truth);
    free(rhok_1truth);
    free(Fk);

    MPI_Finalize();
    return 0;
}

