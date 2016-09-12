#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/constrainedgaussian.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_void_msg_handler, comm, NULL);

    FastPMSolver * solver = & (FastPMSolver) {
        .nc = 128,
        .boxsize = 32.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
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

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 2e6, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    FastPM2PCF xi;

    fastpm_2pcf_from_powerspectrum(&xi, (fastpm_fkfunc) fastpm_utils_powerspec_eh, &eh);

    FastPMConstrainedGaussian cg = {
            .constraints = (FastPMConstraint[]) {
                {{32, 32, 32}, 10},
                {{64, 64, 64}, 10},
                {{-1, -1, -1}, -1},
            },
            .size = 2, /* get rid of this.*/
    };

    fastpm_cg_induce_correlation(&cg, solver->basepm, &xi, rho_init_ktruth);

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

