#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>
#include <fastpm/fof.h>
#include <fastpm/io.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 64,
        .boxsize = 64.,
        .alloc_factor = 10.0,
        .cosmology = NULL,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 1},
            {.a_start = 0.0001, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
    };

    FastPMSolver solver[1];
    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 1e7, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->pm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, NULL, 0.1);

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));


    FastPMFOFFinder fof = {
        .linkinglength = 0.2,
        .nmin = 32,
        .kdtree_thresh = 8,
    };

    fastpm_fof_init(&fof, fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM), solver->pm);

    FastPMStore halos[1];

    fastpm_store_set_name(halos, "FOFHalos");
    fastpm_fof_execute(&fof, halos);

    char * snapshot = fastpm_strdup_printf("fof-%d", solver->NTask);
    fastpm_sort_snapshot(halos, solver->comm, FastPMSnapshotSortByLength, 0);
    fastpm_store_write(halos, snapshot, "w", 1, solver->comm);

    int task;
    int ntask;
    MPI_Comm_size(MPI_COMM_WORLD, &ntask);
    MPI_Comm_rank(MPI_COMM_WORLD, &task);
    int j;
    for(j = 0; j < ntask; j ++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(j != task) continue;
        int i;
        for(i = 0; i < 10 && i < halos->np; i ++) {
            fastpm_ilog(INFO, "Length of halo %d: %d\n", i, halos->length[i]);
        }
    }
    fastpm_store_destroy(halos);
    fastpm_fof_destroy(&fof);

    pm_free(solver->pm, rho_final_xtruth);
    pm_free(solver->pm, rho_final_ktruth);
    pm_free(solver->pm, rho_init_ktruth);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

