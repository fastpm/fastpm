#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_void_msg_handler, comm, NULL);

    FastPM * solver = & (FastPM) {
        .nc = 128,
        .boxsize = 32.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .USE_COLA = 0,
        .USE_NONSTDDA = 0,
        .USE_MODEL = 0,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
    };

    fastpm_init(solver, 0, 0, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm_2lpt);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->pm_2lpt);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->pm_2lpt);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->pm_2lpt, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm_2lpt, rho_init_ktruth, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh);

    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
    fastpm_setup_ic(solver, rho_init_ktruth);

    fastpm_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    fastpm_utils_paint(solver->pm_2lpt, solver->p, rho_final_xtruth, rho_final_ktruth, NULL, 0);
    fastpm_utils_dump(solver->pm_2lpt, "fastpm_rho_final_xtruth.raw", rho_final_xtruth);

    pm_free(solver->pm_2lpt, rho_final_xtruth);
    pm_free(solver->pm_2lpt, rho_final_ktruth);
    pm_free(solver->pm_2lpt, rho_init_ktruth);

    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

