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
        .boxsize = 1600.,
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
    
    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
    
    //ADD NCDM
    double m_ncdm[3] = {0.05, 0.02, 0.01};
    int n_ncdm = 3;
    FastPMncdmInitData* nid = fastpm_ncdm_init_create(config->boxsize, m_ncdm, n_ncdm, 9., 10, 2, 0);
    
    int f_subsample_1d = 4;
    int f_subsample_3d = f_subsample_1d*f_subsample_1d*f_subsample_1d;
    FastPMStore * cdm = fastpm_solver_get_species(solver, FASTPM_SPECIES_CDM);

    size_t total_np_cdm = fastpm_store_get_np_total(cdm, comm);

    size_t total_np_ncdm = total_np_cdm / f_subsample_3d * nid->n_split;

    FastPMStore ncdm[1];
    fastpm_store_init_evenly(ncdm,
          fastpm_species_get_name(FASTPM_SPECIES_NCDM),
          total_np_ncdm,
          COLUMN_POS | COLUMN_VEL | COLUMN_ID | COLUMN_MASK | COLUMN_MASS | COLUMN_ACC | solver->config->ExtraAttributes,
          solver->config->alloc_factor,
          solver->comm);

    fastpm_solver_add_species(solver, 
                              FASTPM_SPECIES_NCDM, 
                              ncdm);

    //END OF ADD NCDM

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


    fastpm_solver_setup_lpt(solver, FASTPM_SPECIES_CDM, rho_init_ktruth, 0.1);
    

    //SPLIT
    fastpm_split_ncdm(nid, cdm, ncdm, comm);

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

    fastpm_store_destroy(ncdm);
    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

