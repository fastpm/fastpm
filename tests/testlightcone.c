#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/lightcone.h>

int
write_snapshot(FastPMSolver * fastpm, FastPMStore * p, const char * filebase, char * parameters, int Nwriters);

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 128,
        .boxsize = 128.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
        .COMPUTE_POTENTIAL = 1,
    };
    FastPMSolver solver[1];
    FastPMDriftFactor drift;
    FastPMKickFactor kick;
    FastPMLightCone lc[1];

    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5e6, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_info("dx1  : %g %g %g %g\n", 
            solver->info.dx1[0], solver->info.dx1[1], solver->info.dx1[2],
            (solver->info.dx1[0] + solver->info.dx1[1] + solver->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n", 
            solver->info.dx2[0], solver->info.dx2[1], solver->info.dx2[2],
            (solver->info.dx2[0] + solver->info.dx2[1] + solver->info.dx2[2]) / 3.0);
    double time_step[] = {0.1};
    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    fastpm_lc_init(lc, 0.02, solver->cosmology, solver->p);

    double a, d;
    for(a = 0.1; a < 1.0; a += 0.1) {
        d = fastpm_lc_horizon(lc, a);
        fastpm_info("a = %0.04f z = %0.08f d = %g\n", a, 1 / a - 1, d);
    }

    fastpm_drift_init(&drift, solver, 0.1, 0.1, 1.0);
    fastpm_kick_init(&kick, solver, 0.1, 0.1, 1.0);

    fastpm_lc_intersect(lc, &drift, &kick, solver->p);
    fastpm_info("%td particles are in the light cone\n", lc->p->np);

    write_snapshot(solver, lc->p, "lightconeresult", "", 1);

    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    double time_step2[] = {0.1, 1.0};
    fastpm_solver_evolve(solver, time_step2, sizeof(time_step2) / sizeof(time_step2[0]));
    write_snapshot(solver, solver->p, "nonlightconeresultZ=0", "", 1);
    
    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    double time_step3[] = {0.1};
    fastpm_solver_evolve(solver, time_step3, sizeof(time_step3) / sizeof(time_step3[0]));
    write_snapshot(solver, solver->p, "nonlightconeresultZ=9", "", 1);

    fastpm_lc_destroy(lc);

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

