#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
#include <fastpm/lightcone.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

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

    FastPMDrift drift;
    fastpm_drift_init(&drift, solver, 0.1, 0.2, 0.3);

    FastPMLightCone lc[1];

    Cosmology CP = {
        .OmegaM = 0.3,
        .OmegaLambda = 0.7,
    };

    fastpm_lc_init(lc, CP, 10000);

    double a, d;
    for(a = 0.1; a < 1.0; a += 0.1) {
        d = fastpm_lc_horizon(lc, a);
        fastpm_info("a = %0.04f z = %0.08f d = %g\n", a, 1 / a - 1, d * HubbleDistance);
    }

    // fastpm_lc_intersect test
    double solution;
    int status;

    status = fastpm_lc_intersect(lc, &solution, 0.2, 0.5);
    status = fastpm_lc_intersect(lc, &solution, 0.2, -0.5);

    status = fastpm_lc_intersect(lc, &solution, 0.3, -0.1);
    status = fastpm_lc_intersect(lc, &solution, 0.1, -0.3);
    status = fastpm_lc_intersect(lc, &solution, 0.3, 0.1);
    status = fastpm_lc_intersect(lc, &solution, 0.1, 0.3);

    status = fastpm_lc_intersect(lc, &solution, 0.1, 0.1);
    status = fastpm_lc_intersect(lc, &solution, 0.0, 0.1);
    status = fastpm_lc_intersect(lc, &solution, 0.1, 0.0);
    status = fastpm_lc_intersect(lc, &solution, 0.0, 0.0);

    if(status == 0) {
        fastpm_info("last calculation failed");
    }
    fastpm_lc_destroy(lc);

    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

