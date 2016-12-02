#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <stdlib.h>

#include <mpi.h>
#include <pfft.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>



int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    int ThisTask;
    int NTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);
    FastPMMesh mesh[1];
    fastpm_mesh_init(mesh, 1.0, 4, 3, (ptrdiff_t []) {1, 1}, comm);

    FastPMMeshIter * iter;
    iter = fastpm_mesh_iter(mesh, FASTPM_MESH_ITER_K);
    float ** table = fastpm_mesh_iter_make_table(iter, FastPMMeshK4Point);

    fastpm_info("%g %g \n", table[0][0], table[0][1]);
    int count = 0; 
    count = 0;
    while(fastpm_mesh_iter_next(iter)) {
        fastpm_ilog(0, "%td %td %td %td k4point = %g %g %g\n", iter->iabs[0], iter->iabs[1], iter->iabs[2],
            iter->ind,
            table[0][iter->iabs[0]],
            table[1][iter->iabs[1]],
            table[2][iter->iabs[2]]);
        count++;
    }
    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, comm);
    if(count != 4 * 4 * 3) {
        fastpm_raise(-1, "Total count of the K iterator is wrong\n");
    }
    free(table);
    fastpm_mesh_iter_free(iter);

    iter = fastpm_mesh_iter(mesh, FASTPM_MESH_ITER_X);
    count = 0;
    while(fastpm_mesh_iter_next(iter)) {
        fastpm_ilog(0, "%td %td %td %td\n",
            iter->iabs[0], iter->iabs[1], iter->iabs[2],
            iter->ind);
        count ++;
    }
    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, comm);
    if(count != 4 * 4 * 4) {
        fastpm_raise(-1, "Total count of the K iterator is wrong\n");
    }
    fastpm_mesh_iter_free(iter);

    fastpm_mesh_destroy(mesh);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

