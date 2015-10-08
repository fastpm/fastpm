#include <stdio.h>
#include "libfastpm.h"

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    msg_init(comm);
    msg_set_loglevel(verbose);
    timer_init();

    int nc = 128;
    double boxsize = 256;
    int NTask;
    PMStore pdata;
    PM pm;

    MPI_Comm_size(comm, &NTask);

    pm_store_init(&pdata);

    pm_store_alloc(&pdata, 1.0 * pow(nc, 3) / NTask * 2);

    pm_2lpt_init(&pm, &pdata, nc, boxsize, comm);

    pm_start(&pm);

    /* fill in the modes */
    pm_2lpt_main(&pm, &pdata, comm);

    pm_destroy(&pm);
    MPI_Finalize();
    return 0;
}

