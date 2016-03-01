#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/transfer.h>

static void 
get_lagrangian_position(PMStore * p, ptrdiff_t index, double pos[3]) 
{
    pos[0] = p->q[index][0];
    pos[1] = p->q[index][1];
    pos[2] = p->q[index][2];
}
static void 
get_position(PMStore * p, ptrdiff_t index, double pos[3]) 
{
    pos[0] = p->x[index][0];
    pos[1] = p->x[index][1];
    pos[2] = p->x[index][2];
}

void 
fastpm_za_hmc_force(FastPM2LPTSolver * solver,
        FastPMFloat * data_x, /* rhop in x-space*/
        FastPMFloat * sigma_x, /* sigma_x in x-space*/
        FastPMFloat * Fk,    /* (out) hmc force in fourier space */
        double sml
        )
{
    int d;

    FastPMFloat * workspace = pm_alloc(solver->pm);
    FastPMFloat * workspace2 = pm_alloc(solver->pm);

    /* First obtain a smoothed version of model */
    fastpm_utils_paint(solver->pm, solver->p, NULL, Fk, get_position, 0);
    fastpm_apply_smoothing_transfer(solver->pm, Fk, workspace, sml);
    pm_c2r(solver->pm, workspace);
    ptrdiff_t ind;
    for(ind = 0; ind < pm_size(solver->pm); ind ++) {
        workspace[ind] -= data_x[ind];
        if(sigma_x)
            workspace[ind] /= sigma_x[ind] * sigma_x[ind];
    }

    pm_r2c(solver->pm, workspace, Fk);

    /* Fk contains rhod_k at this point */

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) { 
        fastpm_apply_diff_transfer(solver->pm, Fk, workspace, d);

        /* workspace stores \Gamma(k) = -i k \rho_d */

        pm_c2r(solver->pm, workspace);
        
        fastpm_utils_readout(solver->pm, solver->p, workspace, get_position, ACC[d]);
        fastpm_info("d= %d ACC = %d\n", d, ACC[d]);
    }

    /* now we paint \Psi by the lagrangian position q */

    memset(Fk, 0, sizeof(Fk[0]) * pm_size(solver->pm));

    for(d = 0; d < 3; d ++) {
        fastpm_utils_paint(solver->pm, solver->p, NULL, workspace2, get_lagrangian_position, ACC[d]);
        fastpm_apply_za_hmc_force_transfer(solver->pm, workspace2, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm_size(solver->pm); ind ++) {
            /* Wang's magic factor of 2 in 1301.1348. 
             * We do not put it in in hmc_force_2lpt_transfer */
            Fk[ind] += - 2 * workspace[ind]; 
        }
    }
    pm_free(solver->pm, workspace2);
    pm_free(solver->pm, workspace);
}


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();
    MPI_Comm comm = MPI_COMM_WORLD;
    printf("%td\n", sizeof(FastPMFloat));
    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPM2LPTSolver * solver = & (FastPM2LPTSolver) {
        .USE_DX1_ONLY = 1,
    };


    int nc = 64;
    double boxsize = 512.;
    double alloc_factor = 2.0;
    double omega_m = 0.292;
    int mode = atoi(argv[1]);
    printf("mode = %d %s\n", mode, argv[1]);
    fastpm_2lpt_init(solver, nc, boxsize, alloc_factor, comm);

    FastPMFloat * sigma = pm_alloc(solver->pm);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->pm);
    FastPMFloat * Fk = pm_alloc(solver->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->pm);
    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);

    FastPMFloat * rho_init_k = pm_alloc(solver->pm);
    FastPMFloat * rho_final_k = pm_alloc(solver->pm);
    FastPMFloat * rho_final_x = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5.0e6, /* FIXME: sigma8 ~ 0.8 */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    /*
    double k;
    for(k = 0; k < 1; k += 0.01) {
        double p = fastpm_utils_powerspec_eh(k, &eh);
        printf("%g %g\n", k, p);
    } */

    fastpm_utils_fill_deltak(solver->pm, rho_init_ktruth, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_ktruth[mode * 2] *= 1.02;
    fastpm_2lpt_evolve(solver, rho_init_ktruth, 1.0, omega_m);

    fastpm_utils_paint(solver->pm, solver->p, rho_final_xtruth, rho_final_ktruth, NULL, 0);
    fastpm_utils_dump(solver->pm, "rho_init_ktruth1.raw", rho_init_ktruth);
    fastpm_utils_dump(solver->pm, "rho_final_ktruth1.raw", rho_final_ktruth);
    fastpm_utils_dump(solver->pm, "rho_final_xtruth1.raw", rho_final_xtruth);

    fastpm_utils_fill_deltak(solver->pm, rho_init_ktruth, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    fastpm_2lpt_evolve(solver, rho_init_ktruth, 1.0, omega_m);

    fastpm_utils_paint(solver->pm, solver->p, rho_final_xtruth, rho_final_ktruth, NULL, 0);
    fastpm_utils_dump(solver->pm, "rho_init_ktruth.raw", rho_init_ktruth);
    fastpm_utils_dump(solver->pm, "rho_final_ktruth.raw", rho_final_ktruth);
    fastpm_utils_dump(solver->pm, "rho_final_xtruth.raw", rho_final_xtruth);

    memset(rho_final_x, 0, sizeof(FastPMFloat) * pm_size(solver->pm));
    fastpm_2lpt_hmc_force(solver, rho_final_x, NULL, Fk, 0.);
    fastpm_utils_dump(solver->pm, "Fk.raw", Fk);
    fastpm_za_hmc_force(solver, rho_final_x, NULL, Fk, 0.);
    fastpm_utils_dump(solver->pm, "Fk2.raw", Fk);

    goto byebye;

    /* now a fake MCMC loop. */
    int mcmcstep;

    for(mcmcstep = 0; mcmcstep < 1; mcmcstep++) {
        /* example looping over ks */
        double fac = sqrt(pow(2 * M_PI, 3) / pow(boxsize, 3));
        PMKIter kiter;
        for(pm_kiter_init(solver->pm, &kiter); 
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {

            int d;
            double k[3];
            double kk = 0.0;
            for(d = 0; d < 3; d ++) {
                /* where to get k and k^2 */
                k[d] = kiter.fac[d][kiter.iabs[d]].k_finite;
                kk += kiter.fac[d][kiter.iabs[d]].kk_finite;
            }
         
            k[d] *= 1.0; /* shut up compiler warning */

            double P = fastpm_utils_powerspec_eh(sqrt(kk), &eh);
            rho_init_k[kiter.ind + 0] = fac * sqrt(P) * 1.0; /* real */
            rho_init_k[kiter.ind + 1] = fac * sqrt(P) * 1.0; /* imag */
        }

        /* Evolve to rho_final_k */
        fastpm_2lpt_evolve(solver, rho_init_k, 1.0, omega_m);
        fastpm_utils_paint(solver->pm, solver->p, rho_final_x, rho_final_k, NULL, 0);

        /* calculate the HMC force into Fk */
        fastpm_2lpt_hmc_force(solver, rho_final_xtruth, NULL, Fk, 0);
    }
byebye:

    pm_free(solver->pm, rho_final_x);
    pm_free(solver->pm, rho_final_k);
    pm_free(solver->pm, rho_init_k);

    pm_free(solver->pm, rho_init_ktruth);
    pm_free(solver->pm, rho_final_xtruth);
    pm_free(solver->pm, Fk);
    pm_free(solver->pm, rho_final_ktruth);
    pm_free(solver->pm, sigma);

    fastpm_2lpt_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

