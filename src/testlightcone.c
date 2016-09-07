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
    fastpm_lc_destroy(lc);

    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

