#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/lc-unstruct.h>
#include <fastpm/io.h>


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 64,
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

    double tiles[1][3] = {
            {0, 0, 0},
            };

    FastPMLightCone lc[1] = {{
        .speedfactor = 0.2,
        .glmatrix = {
                {1, 0, 0, 0,},
                {0, 1, 0, 0,},
                {0, 0, 1, 0,},
                {0, 0, 0, 0,},
            },

        .fov = 0.,
        .cosmology = solver->cosmology,
    }};

    FastPMUSMesh usmesh[1];
    FastPMSMesh  smesh[1];

    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_lc_init(lc);
    fastpm_usmesh_init(usmesh, lc, solver->p->np_upper, tiles, 1);
    {
        double xy[][2] =  {{0, 0}, {1, 1,}};
        double z[] = {0, 1, 2, 3};
        fastpm_smesh_init_plane(smesh, lc, xy, 2, z, 4);
    }
    fastpm_info("dx1  : %g %g %g %g\n",
            solver->info.dx1[0], solver->info.dx1[1], solver->info.dx1[2],
            (solver->info.dx1[0] + solver->info.dx1[1] + solver->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n",
            solver->info.dx2[0], solver->info.dx2[1], solver->info.dx2[2],
            (solver->info.dx2[0] + solver->info.dx2[1] + solver->info.dx2[2]) / 3.0);

    double time_step[] = {0.1};
    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    double a, d;
    for(a = 0.1; a < 1.0; a += 0.1) {
        d = HorizonDistance(a, lc->horizon);
        fastpm_info("a = %0.04f z = %0.08f d = %g\n", a, 1 / a - 1, d);
    }

    fastpm_drift_init(&drift, solver, 0.1, 0.1, 1.0);
    fastpm_kick_init(&kick, solver, 0.1, 0.1, 1.0);

    fastpm_usmesh_intersect(usmesh, &drift, &kick, solver);
    fastpm_info("%td particles are in the light cone\n", usmesh->p->np);

    write_snapshot(solver, usmesh->p, "lightconeresult-p", "", 1, NULL);

    fastpm_smesh_compute_potential(smesh, solver->basepm, solver->gravity, rho_init_ktruth, 1.0);

    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    double time_step2[] = {0.1, 1.0};
    fastpm_solver_evolve(solver, time_step2, sizeof(time_step2) / sizeof(time_step2[0]));
    write_snapshot(solver, solver->p, "nonlightconeresultZ=0", "", 1, NULL);

    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    double time_step3[] = {0.1};
    fastpm_solver_evolve(solver, time_step3, sizeof(time_step3) / sizeof(time_step3[0]));
    write_snapshot(solver, solver->p, "nonlightconeresultZ=9", "", 1, NULL);

    fastpm_smesh_destroy(smesh);
    fastpm_usmesh_destroy(usmesh);
    fastpm_lc_destroy(lc);

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_solver_destroy(solver);

    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
