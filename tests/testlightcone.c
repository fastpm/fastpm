#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/io.h>
#include <fastpm/string.h>

static void
interp_handler(FastPMSolver * fastpm, FastPMInterpolationEvent * event, FastPMUSMesh * usmesh)
{
    fastpm_usmesh_intersect(usmesh, event->drift, event->kick,
            event->drift->ai,
            event->drift->af,
            event->whence, fastpm->comm);
}

double tiles[4*4*4][3];
double a[128];

static void
stage1(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{

    FastPMUSMesh usmesh[1];

    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, NULL, 0.1);

    FastPMStore * p = fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM);
    fastpm_usmesh_init(usmesh, lc, pm_volume(solver->pm), p, p->np_upper, tiles, sizeof(tiles) / sizeof(tiles[0]), 0.4, 0.8);

    fastpm_info("stage 1\n");

    double time_step[] = {0.1};

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    FastPMDriftFactor drift;
    FastPMKickFactor kick;

    fastpm_drift_init(&drift, solver, 0.1, 0.1, 1.0);
    fastpm_kick_init(&kick, solver, 0.1, 0.1, 1.0);

    fastpm_usmesh_intersect(usmesh, &drift, &kick, 0.1, 1.0, TIMESTEP_CUR, solver->comm);
    fastpm_info("%td particles are in the light cone\n", usmesh->p->np);

    fastpm_store_write(usmesh->p, "lightconeresult-p", "w", 1, solver->comm);

    FastPMPainter painter[1];

    fastpm_painter_init(painter, solver->pm, solver->config->PAINTER_TYPE, solver->config->painter_support);
    fastpm_usmesh_destroy(usmesh);

}

static void
stage2(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{
    fastpm_info("stage 2\n");
    FastPMUSMesh usmesh[1];

    FastPMStore * p = fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM);
    fastpm_usmesh_init(usmesh, lc, pm_volume(solver->pm), p, p->np_upper, tiles, sizeof(tiles) / sizeof(tiles[0]), 0.4, 0.8);

    fastpm_add_event_handler(&solver->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) interp_handler,
        usmesh);


    double time_step2[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, NULL, 0.1);

    fastpm_solver_evolve(solver, time_step2, sizeof(time_step2) / sizeof(time_step2[0]));

    fastpm_store_write(usmesh->p, "lightcone-unstruct", "w", 1, solver->comm);

    fastpm_remove_event_handler(&solver->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) interp_handler,
        usmesh);

    fastpm_usmesh_destroy(usmesh);
}

static void
stage3(FastPMSolver * solver, FastPMLightCone * lc, FastPMFloat * rho_init_ktruth)
{
    fastpm_info("stage 3\n");

    FastPMUSMesh usmesh[1];

    FastPMStore * p = fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM);
    fastpm_usmesh_init(usmesh, lc, pm_volume(solver->pm), p, p->np_upper, tiles, sizeof(tiles) / sizeof(tiles[0]), 0.4, 0.8);

    fastpm_add_event_handler(&solver->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) interp_handler,
        usmesh);

    double time_step3[] = {0.1};
    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, NULL, 0.1);
    fastpm_solver_evolve(solver, time_step3, sizeof(time_step3) / sizeof(time_step3[0]));

    fastpm_store_write(fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM),
            "nonlightconeresultZ=9", "w", 1, solver->comm);

    fastpm_remove_event_handler(&solver->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) interp_handler,
        usmesh);

    fastpm_usmesh_destroy(usmesh);
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
        .cosmology = NULL,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .ExtraAttributes = COLUMN_POTENTIAL,
    };
    FastPMSolver solver[1];

    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5e6, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->pm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

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

    pm_free(solver->pm, rho_init_ktruth);
    fastpm_solver_destroy(solver);

    libfastpm_cleanup();
    MPI_Finalize();

    return 0;
}
