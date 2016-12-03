#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
    fastpm_column_init_float(x, 3, 128);

    FastPMColumn h[1];
    fastpm_column_init_double_const(h, 3, (double[]){0.5, 0.5, 0.5});

    FastPMColumn f[1];
    fastpm_column_init_float(f, 3, 128);

    FastPMMesh mesh[1];

    FastPMDomain domain[1];
    fastpm_mesh_init(mesh, 3, 4., 4, NULL, comm);

    size_t oldsize = 4;
    fastpm_column_resize(x, oldsize);
    fastpm_column_resize(h, oldsize);
    fastpm_column_resize(f, oldsize);
    int i;

    for(i = 0; i < oldsize; i ++) {
        double pos[3] = {i % 4, i % 4, i % 4};
        fastpm_column_set_double(x, i, pos);
        fastpm_column_set_double(f, i, pos);
    }

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }

    for(i = 0; i < x->size; i ++) {
        double pos[3];
        int ipos[3];
        fastpm_column_get_double(x, i, pos);
        int master;
        int d;
        for(d = 0; d < 3; d ++)
            ipos[d] = fastpm_mesh_pos_to_ipos(mesh, pos[d], d);
        master = fastpm_mesh_ipos_to_rank(mesh, ipos);
        printf("ThisTask = %d pos[%d] =%d %d %d master = %d\n", ThisTask, i, ipos[0], ipos[1], ipos[2], master);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_info("testing domain\n");
    fastpm_domain_init(domain, mesh, x, comm);

    fastpm_domain_decompose(domain, (FastPMColumn * []) {x, h, NULL});

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    for(i = 0; i < x->size; i ++) {
        double pos[3];
        int ipos[3];
        fastpm_column_get_double(x, i, pos);
        int master;
        int d;
        for(d = 0; d < 3; d ++)
            ipos[d] = fastpm_mesh_pos_to_ipos(mesh, pos[d], d);
        master = fastpm_mesh_ipos_to_rank(mesh, ipos);
        printf("ThisTask = %d pos[%d] =%g %g %g master = %d\n", ThisTask, i, pos[0], pos[1], pos[2], master);
        if(master != ThisTask) abort();
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }
    fastpm_info("testing ghosts\n");
    FastPMGhosts ghosts[1];

    fastpm_ghosts_init(ghosts, mesh, x, h, comm);

    FastPMColumn * ghost_x = fastpm_ghosts_fetch(ghosts, x);
    FastPMColumn * ghost_h = fastpm_ghosts_fetch(ghosts, h);
    FastPMColumn * ghost_f = fastpm_ghosts_fetch(ghosts, f);

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    for(i = 0; i < NTask; i ++) {
        printf("sendcount[%d] = %d\n", i, ghosts->sendcounts[i]);
        printf("recvcount[%d] = %d\n", i, ghosts->recvcounts[i]);
    }
    for(i = 0; i < ghost_x->size; i ++) {
        double pos[3] = {0, 0, 0};
        int ipos[3];
        double margin[3] = {0, 0, 0};
        fastpm_column_get_double(ghost_x, i, pos);
        fastpm_column_get_double(ghost_h, i, margin);
        int master;
        int d;
        for(d = 0; d < 3; d ++)
            ipos[d] = fastpm_mesh_pos_to_ipos(mesh, pos[d], d);
        master = fastpm_mesh_ipos_to_rank(mesh, ipos);
        printf("ThisTask = %d pos[%d] =%g %g %g  margin = %g master = %d\n", ThisTask, i, pos[0], pos[1], pos[2], margin[0], master);
        if(master == ThisTask) abort();
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }
    fastpm_info("testing ghost reduction\n");

    for(i = 0; i < x->size; i ++) {
        fastpm_column_set_double(f, i, (double[]) {1.0, 1.0, 1.0});
    }

    for(i = 0; i < ghost_x->size; i ++) {
        fastpm_column_set_double(ghost_f, i, (double[]) {1.0, 1.0, 1.0});
    }
    fastpm_ghosts_reduce(ghosts, f, ghost_f);

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("x->size = %td\n", x->size);
    for(i = 0; i < x->size; i ++) {
        double pos[3] = {0, 0, 0};
        double force[3] = {0, 0, 0};
        double margin[3] = {0, 0, 0};
        fastpm_column_get_double(x, i, pos);
        fastpm_column_get_double(h, i, margin);
        fastpm_column_get_double(f, i, force);
        printf("ThisTask = %d pos[%d] =%g %g %g  f = %g \n", ThisTask, i, pos[0], pos[1], pos[2], force[0]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_free(ghost_f);
    fastpm_column_free(ghost_h);
    fastpm_column_free(ghost_x);

    fastpm_ghosts_destroy(ghosts);
    fastpm_domain_destroy(domain);
    fastpm_mesh_destroy(mesh);
    fastpm_column_destroy(f);
    fastpm_column_destroy(h);
    fastpm_column_destroy(x);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
