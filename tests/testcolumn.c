#include <stdio.h>
#include <string.h>
#include <alloca.h>

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

    FastPMColumn x[1];
    FastPMColumn v[1];
    int64_t index[128];
    fastpm_column_init_float(x, 3, 256);
    fastpm_column_init_float(v, 3, 256);

    size_t oldsize = (1 + ThisTask) * 10;
    size_t newsize = (NTask - ThisTask) * 10;
    fastpm_column_resize(x, oldsize);
    fastpm_column_resize(v, oldsize);

    FastPMColumnSet p[1];
    fastpm_columnset_init(p, (FastPMColumn *[]){x, v, NULL});

    int i;

    for(i = 0; i < oldsize; i ++) {
        double pos[3] = {ThisTask * 10 + i, ThisTask * 10 + i, i};
        fastpm_column_set_double(x, i, pos);
        fastpm_column_set_double(v, i, pos);
        index[i] = i;
    }

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("oldsize = %td\n", oldsize);
    for(i = 0; i < oldsize; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("oldpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_parallel_permute((FastPMColumn*)p, index, newsize, comm);

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("newsize = %td\n", newsize);
    for(i = 0; i < newsize; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("newpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
        fastpm_column_get_double(v, i, pos);
        printf("newvel[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_destroy((FastPMColumn*) p);
    fastpm_column_destroy(v);
    fastpm_column_destroy(x);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
