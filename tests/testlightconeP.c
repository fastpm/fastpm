#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/lightcone_potential.h>
#include <fastpm/io.h>


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 128,
        .boxsize = 2048,//128.,
        .alloc_factor = 2.0,
        .Omega_m = 0.292,
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

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5e6, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->pm, rho_init_ktruth, 2005, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    FastPMLightConeP lcp[1] = {{
        .speedfactor = 1,//0.2,
        .compute_potential = 1,
        .glmatrix = {
                {1, 0, 0, 0,},
                {0, 1, 0, 0,},
                {0, 0, 1, 0,},
                {0, 0, 0, 0,},
            },
        .read_ra_dec=2,
        .fov = 360.,
        .cosmology = solver->cosmology,
    }};
    double (*tiles)[3];
    int ntiles;

    if (lcp->read_ra_dec!=0){
      tiles=fastpm_lcp_tile(solver,0, 0, 0, 0, &ntiles, tiles);
      lcp->shell_params=malloc(sizeof(grid_params)*lcp->read_ra_dec);
      lcp->shell_params[0].ra_dec_filename="healpix512";
      lcp->shell_params[1].ra_dec_filename="healpix64";
      lcp->shell_params[0].a_min=.5;
      lcp->shell_params[0].a_max=.9;

      lcp->shell_params[1].a_min=.9;
      lcp->shell_params[1].a_max=1;
      lcp->shell_params[0].n_a=25;
      lcp->shell_params[1].n_a=5;
      lcp->shell_params[0].subsample_factor=1;
      lcp->shell_params[1].subsample_factor=1;
    }
    else{
      tiles=fastpm_lcp_tile(solver, 1, 1, 1, 1, &ntiles, tiles);
    }


    fastpm_solver_setup_lpt(solver, rho_init_ktruth);

    fastpm_lcp_init(lcp, solver, tiles, ntiles);

    fastpm_info("dx1  : %g %g %g %g\n",
            solver->info.dx1[0], solver->info.dx1[1], solver->info.dx1[2],
            (solver->info.dx1[0] + solver->info.dx1[1] + solver->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n",
            solver->info.dx2[0], solver->info.dx2[1], solver->info.dx2[2],
            (solver->info.dx2[0] + solver->info.dx2[1] + solver->info.dx2[2]) / 3.0);
    // double time_step[] = {0.1,1};
    // fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));
    //
    // double a, d;
    // for(a = 0.1; a < 1.0; a += 0.1) {
    //     d = fastpm_lcp_horizon(lcp, a);
    //     fastpm_info("a = %0.04f z = %0.08f d = %g\n", a, 1 / a - 1, d);
    // }
    //
    // fastpm_drift_init(&drift, solver, 0.1, 0.1, 1.0);
    // fastpm_kick_init(&kick, solver, 0.1, 0.1, 1.0);
    //
    // //fastpm_lcp_intersect(lcp, &drift, &kick, solver);
    // fastpm_info("%td particles are in the light cone\n", lcp->p->np);
    // fastpm_info("%td uniform particles are in the light cone\n", lcp->q->np);
    //
    // write_snapshot(solver, lcp->p, "lightconePresult-p", "", 1, NULL);
    // write_snapshot(solver, lcp->q, "lightconePresult-q", "", 1, NULL);

    // fastpm_solver_setup_lpt(solver, rho_init_ktruth);
    // fastpm_lcp_init(lcp, solver, tiles, ntiles);

    double time_step2[] = {0.1,0.5,0.7,0.9,.99,1};
    fastpm_solver_evolve(solver, time_step2, sizeof(time_step2) / sizeof(time_step2[0]));
    fastpm_info("%td particles are in the light cone\n", lcp->p->np);
    fastpm_info("%td uniform particles are in the light cone\n", lcp->q->np);

    write_snapshot(solver, solver->p, "nonlightconePresultZ=0", "", 1, NULL);

    write_snapshot(solver, lcp->p, "lightconePresult-p2", "", 1, NULL);
    write_snapshot(solver, lcp->q, "lightconePresult-q2", "", 1, NULL);

    // fastpm_solver_setup_lpt(solver, rho_init_ktruth);
    // fastpm_lcp_init(lcp, solver, tiles, ntiles);
    // double time_step3[] = {0.1};
    // fastpm_solver_evolve(solver, time_step3, sizeof(time_step3) / sizeof(time_step3[0]));
    // write_snapshot(solver, solver->p, "nonlightconeresultZ=9", "", 1, NULL);

    fastpm_lcp_destroy(lcp);

    pm_free(solver->pm, rho_init_ktruth);
    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
