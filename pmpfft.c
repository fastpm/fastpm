#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "pmpfft.h"

static MPI_Datatype MPI_PTRDIFF = NULL;


static void module_init() {
    if(MPI_PTRDIFF != NULL) return;
        
    if(sizeof(ptrdiff_t) == 8) {
        MPI_PTRDIFF = MPI_LONG;
    } else {
        MPI_PTRDIFF = MPI_INT;
    }
}

void pm_pfft_init(PM * pm, PMInit * init, PMIFace * iface, MPI_Comm comm) {

    module_init();

    pm->init = *init;
    pm->iface = *iface;

    /* initialize the domain */    
    MPI_Comm_rank(comm, &pm->ThisTask);
    MPI_Comm_size(comm, &pm->NTask);

    int Nx = init->NprocX;
    int Ny;
    if(Nx <= 0) {
        Nx = 1;
        Ny = pm->NTask;
        for(; Nx * Nx < pm->NTask; Nx ++) continue;
        for(; Nx >= 1; Nx--) {
            if (pm->NTask % Nx == 0) break;
            continue;
        }
    } else {
        if(pm->NTask % Nx != 0) {
            fprintf(stderr, "NprocX(%d) and NTask(%d) is incompatible\n", Nx, pm->NTask);
            MPI_Abort(comm, -1);
        }
    }
    Ny = pm->NTask / Nx;
    pm->Nproc[0] = Nx;
    pm->Nproc[1] = Ny;

    pm->Nmesh[0] = init->Nmesh;
    pm->Nmesh[1] = init->Nmesh;
    pm->Nmesh[2] = init->Nmesh;

    pm->BoxSize[0] = init->BoxSize;
    pm->BoxSize[1] = init->BoxSize;
    pm->BoxSize[2] = init->BoxSize;

    pm->Below[0] = 0;
    pm->Below[1] = 0;
    pm->Below[2] = 0;

    pm->Above[0] = 1;
    pm->Above[1] = 1;
    pm->Above[2] = 1;

    pfft_create_procmesh(2, comm, pm->Nproc, &pm->Comm2D);
    pm->allocsize = 2 * pfft_local_size_dft_r2c(
                3, pm->Nmesh, pm->Comm2D, PFFT_TRANSPOSED_OUT | PFFT_PADDED_R2C, 
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);


    /* Set up strides for IRegion (real) and ORegion(complex) */
    pm->IRegion.strides[2] = 1;
    pm->IRegion.strides[1] = 2* (pm->Nmesh[2] / 2 + 1); /* padded */
    pm->IRegion.strides[0] = pm->IRegion.size[1] * pm->IRegion.strides[1];

    pm->ORegion.strides[2] = 1;
    pm->ORegion.strides[0] = pm->Nmesh[2] / 2 + 1; /* transposed */
    pm->ORegion.strides[1] = pm->IRegion.size[0] * pm->IRegion.strides[0];

    int d;
    for(d = 0; d < 2; d ++) {
        MPI_Comm projected;
        int remain_dims[2] = {0, 0};
        remain_dims[d] = 1; 

        pm->Grid.edges_int[d] = 
            malloc(sizeof(pm->Grid.edges_int[0][0]) * (pm->Nproc[d] + 1));
        pm->Grid.edges_float[d] = 
            malloc(sizeof(pm->Grid.edges_float[0][0]) * (pm->Nproc[d] + 1));

        pm->Grid.MeshtoCart[d] = malloc(sizeof(int) * pm->Nmesh[d]);

        MPI_Cart_sub(pm->Comm2D, remain_dims, &projected);
        MPI_Allgather(&pm->IRegion.start[d], 1, MPI_PTRDIFF, 
            pm->Grid.edges_int[d], 1, MPI_PTRDIFF, projected);
        int ntask;
        MPI_Comm_size(projected, &ntask);

        MPI_Comm_free(&projected);
        int j;
        for(j = 0; j < pm->Nproc[d]; j ++) {
            pm->Grid.edges_float[d][j] = 1.0 * pm->Grid.edges_int[d][j] / pm->Nmesh[d] * pm->BoxSize[d];
        }
        /* Last edge is at the edge of the box */
        pm->Grid.edges_float[d][j] = pm->BoxSize[d];
        pm->Grid.edges_int[d][j] = pm->Nmesh[d];
        /* fill in the look up table */
        for(j = 0; j < pm->Nproc[d]; j ++) {
            int i;
            for(i = pm->Grid.edges_int[d][j]; i < pm->Grid.edges_int[d][j+1]; i ++) {
                pm->Grid.MeshtoCart[d][i] = j;
            }
        }
    }

    void * buffer = pm->iface.malloc(pm->allocsize * sizeof(double));

    pm->r2c = pfft_plan_dft_r2c(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_FORWARD, PFFT_PADDED_R2C | PFFT_TRANSPOSED_OUT | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);

    pm->c2r = pfft_plan_dft_c2r(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_BACKWARD, PFFT_PADDED_C2R | PFFT_TRANSPOSED_IN | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);
    pm->iface.free(buffer);

    for(d = 0; d < 3; d++) {
        pm->MeshtoK[d] = malloc(pm->Nmesh[d] * sizeof(double));
        int i;
        for(i = 0; i < pm->Nmesh[d]; i++) {
            int ii = i;
            if(ii > pm->Nmesh[d] / 2) {
                ii -= pm->Nmesh[d];
            }
            pm->MeshtoK[d][i] = i * 2 * M_PI / pm->BoxSize[d];
        }
    }
}

int pm_pos_to_rank(PM * pm, double pos[3]) {
    int d;
    int rank2d[2];
    for(d = 0; d < 2; d ++) {
        int ipos = floor(pos[d] / pm->BoxSize[d] * pm->Nmesh[d]);
        while(ipos < 0) ipos += pm->Nmesh[d];
        while(ipos >= pm->Nmesh[d]) ipos -= pm->Nmesh[d];
        rank2d[d] = pm->Grid.MeshtoCart[d][ipos];
    }
    return rank2d[0] * pm->Nproc[1] + rank2d[1];
}
void pm_start(PM * pm) {
    pm->canvas = pm->iface.malloc(sizeof(double) * pm->allocsize);
    pm->workspace = pm->iface.malloc(sizeof(double) * pm->allocsize);
}
void pm_stop(PM * pm) {
    pm->iface.free(pm->canvas);
    pm->iface.free(pm->workspace);
    pm->canvas = NULL;
    pm->workspace = NULL;
}

void pm_r2c(PM * pm) {
    pfft_execute_dft_r2c(pm->r2c, pm->canvas, (pfft_complex*)pm->canvas);
}

void pm_c2r(PM * pm) {
    pfft_execute_dft_c2r(pm->c2r, (pfft_complex*) pm->workspace, pm->workspace);
}

