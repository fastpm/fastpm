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
#include <fastpm/string.h>

static void
smesh_handler(FastPMSMesh * mesh, FastPMLCEvent * lcevent, FastPMSolver * solver)
{
    fastpm_info("lc event : a0 = %g a1 = %g, n = %td \n", lcevent->a0, lcevent->a1, lcevent->p->np);
    char * fn = fastpm_strdup_printf("lightcone_struct");
    if(lcevent->is_first) {
        fastpm_info("First iteration\n");
        write_snapshot(solver, lcevent->p, fn, "", 1, FastPMSnapshotSortByAEmit);
    } else {
        fastpm_info("not first iteration\n");
        append_snapshot(solver, lcevent->p, fn, "", 1, FastPMSnapshotSortByAEmit);
    }
    free(fn);
}

static void
force_handler(FastPMSolver * solver, FastPMForceEvent * event, FastPMSMesh * smesh)
{
    fastpm_info("force handler, %g %g\n", event->a_f, event->a_n);
    fastpm_smesh_compute_potential(smesh, event->pm, event->gravity, event->delta_k, event->a_f, event->a_n);
}

static void
interp_handler(FastPMSolver * fastpm, FastPMInterpolationEvent * event, FastPMUSMesh * usmesh)
{
    fastpm_usmesh_intersect(usmesh, event->drift, event->kick, fastpm);
}

double tiles[4*4*4][3];
double a[128];

static void
stage1(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{

    FastPMUSMesh usmesh[1];
    FastPMSMesh  smesh[1];

    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_usmesh_init(usmesh, lc, solver->p->np_upper, tiles, sizeof(tiles) / sizeof(tiles[0]), 0.4, 0.8);

    fastpm_smesh_init(smesh, lc, solver->p->np_upper);
    fastpm_smesh_add_layer_healpix(smesh, 32, a, 64, solver->comm);
    fastpm_smesh_add_layer_healpix(smesh, 16, a + 64, 64, solver->comm);

    fastpm_add_event_handler(&smesh->event_handlers,
            FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEventHandlerFunction) smesh_handler,
            solver);

    fastpm_info("dx1  : %g %g %g %g\n",
            solver->info.dx1[0], solver->info.dx1[1], solver->info.dx1[2],
            (solver->info.dx1[0] + solver->info.dx1[1] + solver->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n",
            solver->info.dx2[0], solver->info.dx2[1], solver->info.dx2[2],
            (solver->info.dx2[0] + solver->info.dx2[1] + solver->info.dx2[2]) / 3.0);

    fastpm_info("stage 1\n");

    double time_step[] = {0.1};

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    FastPMDriftFactor drift;
    FastPMKickFactor kick;

    fastpm_drift_init(&drift, solver, 0.1, 0.1, 1.0);
    fastpm_kick_init(&kick, solver, 0.1, 0.1, 1.0);

    fastpm_usmesh_intersect(usmesh, &drift, &kick, solver);
    fastpm_info("%td particles are in the light cone\n", usmesh->p->np);

    write_snapshot(solver, usmesh->p, "lightconeresult-p", "", 1, NULL);

    fastpm_smesh_compute_potential(smesh, solver->basepm, solver->gravity, rho_init_ktruth, 0.1, 0.5);
    fastpm_smesh_compute_potential(smesh, solver->basepm, solver->gravity, rho_init_ktruth, 0.5, 1.0);
    fastpm_smesh_compute_potential(smesh, solver->basepm, solver->gravity, rho_init_ktruth, 1.0, -1.0);

    fastpm_smesh_destroy(smesh);
    fastpm_usmesh_destroy(usmesh);

}

static void
stage2(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{
    fastpm_info("stage 2\n");
    FastPMUSMesh usmesh[1];
    FastPMSMesh  smesh[1];

    fastpm_usmesh_init(usmesh, lc, solver->p->np_upper, tiles, sizeof(tiles) / sizeof(tiles[0]), 0.4, 0.8);

    fastpm_smesh_init(smesh, lc, solver->p->np_upper);
    fastpm_smesh_add_layers_healpix(smesh, 
            pow(solver->config->nc / solver->config->boxsize, 2),
            pow(solver->config->nc / solver->config->boxsize, 3),
            0.4, 0.8,
            solver->comm);
//    fastpm_smesh_add_layer_healpix(smesh, 32, a, 64, solver->comm);
//    fastpm_smesh_add_layer_healpix(smesh, 16, a + 64, 64, solver->comm);

    fastpm_add_event_handler(&smesh->event_handlers,
            FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEventHandlerFunction) smesh_handler,
            solver);

    fastpm_add_event_handler(&solver->event_handlers,
            FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEventHandlerFunction) force_handler,
            smesh);

    fastpm_add_event_handler(&solver->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) interp_handler,
        usmesh);


    double time_step2[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_solver_evolve(solver, time_step2, sizeof(time_step2) / sizeof(time_step2[0]));

    //write_snapshot(solver, solver->p, "nonlightconeresultZ=0", "", 1, NULL);
    write_snapshot(solver, usmesh->p, "lightcone-unstruct", "", 1, NULL);

    fastpm_remove_event_handler(&solver->event_handlers,
            FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEventHandlerFunction) force_handler,
            smesh);

    fastpm_remove_event_handler(&solver->event_handlers,
            FASTPM_EVENT_INTERPOLATION, FASTPM_EVENT_STAGE_BEFORE,
            (FastPMEventHandlerFunction) interp_handler,
            smesh);

    fastpm_smesh_destroy(smesh);
    fastpm_usmesh_destroy(usmesh);
}

static void
stage3(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{
    fastpm_info("stage 3\n");

    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    double time_step3[] = {0.1};
    fastpm_solver_evolve(solver, time_step3, sizeof(time_step3) / sizeof(time_step3[0]));

    write_snapshot(solver, solver->p, "nonlightconeresultZ=9", "", 1, NULL);

}


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 64,
        .boxsize = 128.,
        .alloc_factor = 10.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .COMPUTE_POTENTIAL = 1,
    };
    FastPMSolver solver[1];

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

    {
        int p = 0;
        int i, j, k;
        for(i = -2; i <= 1; i ++) {
        for(j = -2; j <= 1; j ++) {
        for(k = -2; k <= 1; k ++) {
            tiles[p][0] = i * config->boxsize;
            tiles[p][1] = j * config->boxsize;
            tiles[p][2] = k * config->boxsize;
            p ++;
        }}}
    }

    // for(size_t i = 0; i < npix; i ++) {
    //     fastpm_info("test lightcone %ld %g %g \n",i,ra[i],dec[i]);
    // }

    FastPMLightCone lc[1] = {{
        .speedfactor = 0.01,
        .glmatrix = {
                {0, 1, 0, 0,},
                {1, 0, 0, 0,},
                {0, 0, 1, 0,},
                {0, 0, 0, 1,},
            },

        .fov = 360., /* full sky */
        .cosmology = solver->cosmology,
    }};

    {
        int i;
        for(i = 0; i < 128; i ++) {
            a[i] = 0.4 + 0.4 * i / 128.;
        }
    }

    fastpm_lc_init(lc);

    {
        double a, d;
        for(a = 0.1; a < 1.0; a += 0.1) {
            d = HorizonDistance(a, lc->horizon);
            fastpm_info("a = %0.04f z = %0.08f d = %g\n", a, 1 / a - 1, d);
        }
    }

    stage1(solver, lc, rho_init_ktruth);

    stage2(solver, lc, rho_init_ktruth);

    stage3(solver, lc, rho_init_ktruth);

    fastpm_lc_destroy(lc);

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_solver_destroy(solver);

    libfastpm_cleanup();
    MPI_Finalize();

    return 0;
}
