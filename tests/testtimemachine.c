#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
#include <fastpm/timemachine.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_void_msg_handler, comm, NULL);

    FastPMSolver * solver = & (FastPMSolver) {
        .nc = 128,
        .boxsize = 32.,
        .alloc_factor = 2.0,
        .Omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .USE_NONSTDDA = 0,
        .USE_MODEL = 0,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
    };

    fastpm_init(solver, 0, 0, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    
    fastpm_ic_fill_gaussiank(solver->pm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    double time_step[] = {0.2, 0.4, 0.6, 0.8, 1.0};
    //double time_step[] = {0.1, 1.0};
    fastpm_setup_ic(solver, rho_init_ktruth);

    fastpm_tevo_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->pm, solver->PAINTER_TYPE, solver->painter_support);

    pm_free(solver->pm, rho_init_ktruth);

	fastpm_destroy(solver);

    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

//~ int main(int argc, char * argv[]) {

    //~ MPI_Init(&argc, &argv);

    //~ libfastpm_init();

    //~ MPI_Comm comm = MPI_COMM_WORLD;

    //~ /* timemaschine test */

    //~ entry template[] = {
    //~ {0, 0, 1}, /* Kick */
    //~ {0, 1, 1}, /* Drift */
    //~ {0, 2, 1}, /* Drift */
    //~ {2, 2, 1}, /* Force */
    //~ {2, 2, 2}, /* Kick */
    //~ {-1, -1, -1} /* End of table */
    //~ };
    //~ entry template1[] = {
    //~ {0, 1, 0}, /* Drift */
    //~ {0, 1, 1}, /* Kick */
    //~ {1, 1, 1}, /* Force */
    //~ {-1, -1, -1} /* End of table */
    //~ };

    //~ entry template[] = {
    //~ {0, 0, 1}, /* Kick*/
    //~ {0, 1, 1}, /* Drift*/
    //~ {1, 1, 1}, /* Force */
    //~ {-1, -1, -1} /* End of table */
    //~ };
    
    //~ FastPMTEEntry template[] = {
    //~ {0, 0, 1}, /* Kick */
    //~ {0, 1, 1}, /* Drift */
    //~ {0, 2, 1}, /* Drift */
    //~ {2, 2, 1}, /* Force */
    //~ {2, 2, 2}, /* Kick */
    //~ {-1, -1, -1} /* End of table */
    //~ };

    //~ double timesteps[5] = {0.2, 0.4, 0.6, 0.8, 1.0};
    //~ FastPMTEStates *states = malloc(sizeof(FastPMTEStates));
    
    //~ fastpm_tevo_generate_states(states, 5, template, timesteps);

    //~ fastpm_tevo_print_states(states);
    //~ fastpm_tevo_destroy_states(states);

    //~ /* timemaschine test */

    //~ fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    //~ libfastpm_cleanup();

    //~ MPI_Finalize();
    
    //~ return 0;
//~ }
