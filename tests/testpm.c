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

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 32,
        .boxsize = 32.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
    };

    FastPMSolver solver[1];
    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->basepm);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, 0.1);
    
    //SPLIT
    FastPMncdmInitData* nid = fastpm_ncdm_init_create(0.05, 9., 1, 1);
    //subsample 1/64 = 1/4^3... 4 per dir. first need to build a mask... which partc to keep or not. mask is 8 bit integer, can compute from id... like %2. Similiar routine in store.c create_mask, maybe can reuse here.
    
    int f_subsample = 64;
    size_t np_cdm = fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM)->np;
    size_t np_ncdm = np_cdm / f_subsample * nid->n_split;
    //ADD
    fastpm_solver_add_species(solver, 
                              FASTPM_SPECIES_NCDM, 
                              np_ncdm);
    
    fastpm_split_ncdm(nid, 
                      fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM), 
                      fastpm_solver_get_species(solver, FASTPM_SPECIES_NCDM), 
                      f_subsample);
    fastpm_ncdm_init_free(nid);
    //END SPLIT
    
    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->basepm, config->PAINTER_TYPE, config->painter_support);

    pm_clear(solver->basepm, rho_final_xtruth);
    fastpm_paint(painter, rho_final_xtruth, fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM), FASTPM_FIELD_DESCR_NONE);
    //fastpm_utils_dump(solver->basepm, "fastpm_rho_final_xtruth.raw", rho_final_xtruth);

    pm_free(solver->basepm, rho_final_xtruth);
    pm_free(solver->basepm, rho_final_ktruth);
    pm_free(solver->basepm, rho_init_ktruth);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

