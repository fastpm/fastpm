#include <stdio.h>
#include <string.h>

#include "libfastpm.h"

static void DUMP(FastPM2LPTSolver * solver, char * filename, float_t *data) {
    FILE * fp = fopen(filename, "w");
    fwrite(data, sizeof(float_t), solver->pm->allocsize, fp);
    fclose(fp);
}

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;
    FastPM2LPTSolver * solver = alloca(sizeof(FastPM2LPTSolver));

    int nc = 128;
    double boxsize = 409600.;
    double alloc_factor = 2.0;
    double omega_m = 0.292;

    fastpm_2lpt_init(solver, nc, boxsize, alloc_factor, comm);

    float_t * rho_final_ktruth = pm_alloc(solver->pm);
    float_t * Fk = pm_alloc(solver->pm);
    float_t * rho_final_xtruth = pm_alloc(solver->pm);
    float_t * rho_init_ktruth = pm_alloc(solver->pm);

    float_t * rho_init_k = pm_alloc(solver->pm);
    float_t * rho_final_k = pm_alloc(solver->pm);
    float_t * rho_final_x = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_fill_deltak(solver->pm, rho_init_ktruth, 100, (fastpm_pkfunc)fastpm_powerspec_eh, &eh);

    fastpm_2lpt_evolve(solver, rho_init_ktruth, 1.0, omega_m);

    fastpm_2lpt_paint(solver, rho_final_xtruth, rho_final_ktruth);
    DUMP(solver, "rho_init_ktruth.raw", rho_init_ktruth);
    DUMP(solver, "rho_final_ktruth.raw", rho_final_ktruth);

    memset(rho_final_x, 0, sizeof(float_t) * solver->pm->allocsize);
    fastpm_2lpt_hmc_force(solver, rho_final_x, Fk);
    DUMP(solver, "Fk.raw", Fk);

    /* now a fake MCMC loop. */
    int mcmcstep;

    for(mcmcstep = 0; mcmcstep < 1; mcmcstep++) {
        /* example looping over ks */
        ptrdiff_t ind;
        ptrdiff_t i[3] = {0};
        double fac = sqrt(pow(2 * M_PI, 3) / solver->pm->Volume);

        for(ind = 0; ind < solver->pm->allocsize; ind += 2, pm_inc_o_index(solver->pm, i)) {
            int d;
            double k[3];
            double kk = 0.0;
            for(d = 0; d < 3; d ++) {
                k[d] = solver->pm->MeshtoK[d][i[d] + solver->pm->ORegion.start[d]];
                kk += k[d] * k[d];
            }
         
            double P = fastpm_powerspec_eh(sqrt(kk), &eh);
            rho_init_k[ind + 0] = fac * sqrt(P) * 1.0; /* real */
            rho_init_k[ind + 1] = fac * sqrt(P) * 1.0; /* imag */
        }

        /* Evolve to rho_final_k */
        fastpm_2lpt_evolve(solver, rho_init_k, 1.0, omega_m);
        fastpm_2lpt_paint(solver, rho_final_x, rho_final_k);

        /* calculate the HMC force into Fk */
        fastpm_2lpt_hmc_force(solver, rho_final_xtruth, Fk);
    }

    pm_free(solver->pm, rho_final_x);
    pm_free(solver->pm, rho_final_k);
    pm_free(solver->pm, rho_init_k);

    pm_free(solver->pm, rho_init_ktruth);
    pm_free(solver->pm, rho_final_xtruth);
    pm_free(solver->pm, Fk);
    MPI_Finalize();
    return 0;
}

