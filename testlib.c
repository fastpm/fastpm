#include <stdio.h>
#include <string.h>

#include "libfastpm.h"

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    struct ClockTable ct = {};
    walltime_init(&ct);

    MPI_Comm comm = MPI_COMM_WORLD;
    int nc = 128;
    double boxsize = 256;
    double alloc_factor = 2.0;
    double omega_m = 0.292;

    PMStore * pdata = alloca(sizeof(PMStore));
    PM * pm = alloca(sizeof(PM));

    fastpm_init(pdata, nc, alloc_factor, comm);

    fastpm_init_pm(pm, pdata, nc, boxsize, comm);

    pm_start(pm);

/*
    VPM * vpm_list = vpm_create(3, (int[]) {1, 2, 3}, (double[]){0, 0.2, 0.5},
            &pm->init, &pm->iface, comm);
*/
    float_t * rho_init_k = malloc(sizeof(rho_init_k[0]) * pm->allocsize);
    float_t * rho_final_k = malloc(sizeof(rho_final_k[0]) * pm->allocsize);
    float_t * rho_p_x = malloc(sizeof(rho_p_x[0]) * pm->allocsize);
    float_t * rho_model_x = malloc(sizeof(rho_model_x[0]) * pm->allocsize);
    float_t * rho_init_ktruth = malloc(sizeof(rho_init_ktruth[0]) * pm->allocsize);

    float_t * Fk = malloc(sizeof(Fk[0]) * pm->allocsize);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_fill_deltak(pm, rho_init_ktruth, 100, (fastpm_pkfunc)fastpm_powerspec_eh, &eh);

    fastpm_evolve_2lpt(pm, pdata, 1.0, omega_m, rho_init_ktruth);

/*
    fastpm_evolve_pm(pm, vpm_list, pdata, 
        (double[]) {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 
                10, omega_m, 
            rho_init_ktruth, rho_final_kpm, comm);
*/

    memcpy(rho_p_x, pm->workspace, sizeof(pm->workspace[0]) * pm->allocsize);

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
            rho_init_k[ind + 0] = fac * sqrt(P) * 1.0; /* real */
            rho_init_k[ind + 1] = fac * sqrt(P) * 1.0; /* imag */
        }

        /* Evolve to rho_final_k */
        fastpm_evolve_2lpt(pm, pdata, 1.0, omega_m, rho_init_k);
        memcpy(rho_model_x, pm->workspace, sizeof(pm->workspace[0]) * pm->allocsize);

        pm_r2c(pm);
        memcpy(rho_final_k, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);

        /* calculate the HMC force into Fk */
        fastpm_derivative_2lpt(pm, pdata, rho_p_x, Fk);
    }

    free(rho_init_k);
    free(rho_final_k);
    free(rho_init_ktruth);
    free(Fk);

    pm_stop(pm);
    MPI_Finalize();
    return 0;
}

