#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "pmpfft.h"
#include "msg.h"

static MPI_Datatype MPI_PTRDIFF = NULL;

#if FFT_PRECISION == 64
    #define plan_dft_r2c pfft_plan_dft_r2c
    #define plan_dft_c2r pfft_plan_dft_c2r
    #define execute_dft_r2c pfft_execute_dft_r2c
    #define execute_dft_c2r pfft_execute_dft_c2r
    typedef pfft_plan plan;
    #define _pfft_init pfft_init
#elif FFT_PRECISION == 32
    #define plan_dft_r2c pfftf_plan_dft_r2c
    #define plan_dft_c2r pfftf_plan_dft_c2r
    #define execute_dft_r2c pfftf_execute_dft_r2c
    #define execute_dft_c2r pfftf_execute_dft_c2r
    typedef pfftf_plan plan;
    #define _pfft_init pfftf_init
#endif

static void module_init() {
    if(MPI_PTRDIFF != NULL) return;
        
    _pfft_init();

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
            msg_abort(-1, "NprocX(%d) and NTask(%d) is incompatible\n", Nx, pm->NTask);
        }
    }
    Ny = pm->NTask / Nx;
    pm->Nproc[0] = Nx;
    pm->Nproc[1] = Ny;

    int d;

    pm->Norm = 1.0;
    pm->Volume = 1.0;
    for(d = 0; d < 3; d ++) {
        pm->Nmesh[d] = init->Nmesh;
        pm->BoxSize[d] = init->BoxSize;

        pm->Below[d] = 0;
        pm->Above[d] = 1;

        pm->CellSize[d] = pm->BoxSize[d] / pm->Nmesh[d];
        pm->InvCellSize[d] = 1.0 / pm->CellSize[d]; 
        pm->Norm *= pm->Nmesh[d];
        pm->Volume *= pm->BoxSize[d];
    }


    pfft_create_procmesh(2, comm, pm->Nproc, &pm->Comm2D);
    pm->allocsize = 2 * pfft_local_size_dft_r2c(
                3, pm->Nmesh, pm->Comm2D, 
                          0
                        | (pm->init.transposed?PFFT_TRANSPOSED_OUT:0)
                        | PFFT_PADDED_R2C, 
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);


#if 0
    msg_aprintf(debug, "IRegion : %td %td %td + %td %td %td\n",
        pm->IRegion.start[0],
        pm->IRegion.start[1],
        pm->IRegion.start[2],
        pm->IRegion.size[0],
        pm->IRegion.size[1],
        pm->IRegion.size[2]
    );

    msg_aprintf(debug, "ORegion : %td %td %td + %td %td %td\n",
        pm->ORegion.start[0],
        pm->ORegion.start[1],
        pm->ORegion.start[2],
        pm->ORegion.size[0],
        pm->ORegion.size[1],
        pm->ORegion.size[2]
    );
    /* Set up strides for IRegion (real) and ORegion(complex) */
#endif
    /* Note that we need to fix up the padded size of the real data;
     * and transpose with strides , */


    pm->IRegion.strides[2] = 1;
    pm->IRegion.strides[1] = pm->IRegion.size[2];
    pm->IRegion.strides[0] = pm->IRegion.size[1] * pm->IRegion.strides[1];
    pm->IRegion.total = pm->IRegion.size[0] * pm->IRegion.strides[0];

    pm->IRegion.size[2] = pm->Nmesh[2];

    if(pm->init.transposed) {
        /* transposed, y, z, x */
        pm->ORegion.strides[0] = 1;
        pm->ORegion.strides[2] = pm->ORegion.size[0];
        pm->ORegion.strides[1] = pm->ORegion.size[2] * pm->ORegion.strides[2];
        pm->ORegion.total = pm->ORegion.size[1] * pm->ORegion.strides[1];
    } else {
        /* non-transposed */
        pm->ORegion.strides[2] = 1;
        pm->ORegion.strides[1] = pm->ORegion.size[2];
        pm->ORegion.strides[0] = pm->ORegion.size[1] * pm->ORegion.strides[1];
        pm->ORegion.total = pm->ORegion.size[0] * pm->ORegion.strides[0];
    }

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

    pm->canvas = pm->iface.malloc(pm->allocsize * sizeof(pm->canvas[0]));
    pm->workspace = pm->iface.malloc(pm->allocsize * sizeof(pm->workspace[0]));

    unsigned int r2cflags = PFFT_PADDED_R2C 
            | (pm->init.transposed?PFFT_TRANSPOSED_OUT:0)
            | PFFT_ESTIMATE 
            //| PFFT_MEASURE
            | PFFT_DESTROY_INPUT;
    unsigned int c2rflags = PFFT_PADDED_C2R 
            | (pm->init.transposed?PFFT_TRANSPOSED_IN:0)
            | PFFT_ESTIMATE 
            //| PFFT_MEASURE
            | PFFT_DESTROY_INPUT;

    pm->r2c = (pfft_plan) plan_dft_r2c(
            3, pm->Nmesh, (void*) pm->workspace, (void*) pm->canvas, 
            pm->Comm2D,
            PFFT_FORWARD, r2cflags);

    pm->c2r = (pfft_plan) plan_dft_c2r(
            3, pm->Nmesh, (void*) pm->workspace, (void*) pm->workspace, 
            pm->Comm2D,
            PFFT_BACKWARD, c2rflags);

    pm->iface.free(pm->canvas);
    pm->iface.free(pm->workspace);
    pm->canvas = NULL;
    pm->workspace = NULL;

    for(d = 0; d < 3; d++) {
        pm->MeshtoK[d] = malloc(pm->Nmesh[d] * sizeof(double));
        int i;
        for(i = 0; i < pm->Nmesh[d]; i++) {
            int ii = i;
            if(ii >= pm->Nmesh[d] / 2) {
                ii -= pm->Nmesh[d];
            }
            pm->MeshtoK[d][i] = ii * 2 * M_PI / pm->BoxSize[d];
        }
    }
}

int pm_pos_to_rank(PM * pm, double pos[3]) {
    int d;
    int rank2d[2];
    int ipos[3];
    for(d = 0; d < 2; d ++) {
        ipos[d] = floor(pos[d] * pm->InvCellSize[d]);
    }
    return pm_ipos_to_rank(pm, ipos);
}

int pm_ipos_to_rank(PM * pm, int i[3]) {
    int d;
    int rank2d[2];
    for(d = 0; d < 2; d ++) {
        int ipos = i[d];
        while(UNLIKELY(ipos < 0)) ipos += pm->Nmesh[d];
        while(UNLIKELY(ipos >= pm->Nmesh[d])) ipos -= pm->Nmesh[d];
        rank2d[d] = pm->Grid.MeshtoCart[d][ipos];
    }
    return rank2d[0] * pm->Nproc[1] + rank2d[1];
}

void pm_start(PM * pm) {
    pm->canvas = pm->iface.malloc(sizeof(pm->canvas[0]) * pm->allocsize);
    pm->workspace = pm->iface.malloc(sizeof(pm->canvas[0]) * pm->allocsize);
}
void pm_stop(PM * pm) {
    pm->iface.free(pm->canvas);
    pm->iface.free(pm->workspace);
    pm->canvas = NULL;
    pm->workspace = NULL;
}

void pm_r2c(PM * pm) {
    execute_dft_r2c((plan) pm->r2c, pm->workspace, (void*)pm->canvas);
}

void pm_c2r(PM * pm) {
    execute_dft_c2r((plan) pm->c2r, (void*) pm->workspace, pm->workspace);
}

void pm_unravel_o_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]) {
    /*
     * using pm_unravel_o_index function is slower than pm_inc_o_index, thus it is only used
     * during dev to test pm_inc_o_index.
     * */
    ptrdiff_t tmp = ind;
    if(pm->init.transposed) {
        /* y, z, x*/
        i[1] = tmp / pm->ORegion.strides[1];
        tmp %= pm->ORegion.strides[1];
        i[2] = tmp / pm->ORegion.strides[2];
        tmp %= pm->ORegion.strides[2];
        i[0] = tmp;
    } else {
        i[0] = tmp / pm->ORegion.strides[0];
        tmp %= pm->ORegion.strides[0];
        i[1] = tmp / pm->ORegion.strides[1];
        tmp %= pm->ORegion.strides[1];
        i[2] = tmp;
    }
}
void pm_inc_o_index(PM * pm, ptrdiff_t i[3]) {
    if(pm->init.transposed) {
        i[0] ++;
        if(UNLIKELY(i[0] == pm->ORegion.size[0])) {
            i[0] = 0;
            i[2] ++;
            if(UNLIKELY(i[2] == pm->ORegion.size[2])) {
                i[2] = 0;
                i[1] ++;
            }
        }
    } else {
        i[2] ++;
        if(UNLIKELY(i[2] == pm->ORegion.size[2])) {
            i[2] = 0;
            i[1] ++;
            if(UNLIKELY(i[1] == pm->ORegion.size[1])) {
                i[1] = 0;
                i[0] ++;
            }
        }
    }
}

void pm_unravel_i_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]) {
    ptrdiff_t tmp = ind;
    i[0] = tmp / pm->IRegion.strides[0];
    tmp %= pm->IRegion.strides[0];
    i[1] = tmp / pm->IRegion.strides[1];
    tmp %= pm->IRegion.strides[1];
    i[2] = tmp;
}

void pm_inc_i_index(PM * pm, ptrdiff_t i[3]) {
    i[2] ++;
    if(UNLIKELY(i[2] == pm->IRegion.strides[1])) { /* the padding !*/
        i[2] = 0;
        i[1] ++;
        if(UNLIKELY(i[1] == pm->IRegion.size[1])) {
            i[1] = 0;
            i[0] ++;
        }
    }
}

int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {

    int ThisTask;
    int NTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);
    int PTask;
    int ngrp;

    for(PTask = 0; NTask > (1 << PTask); PTask++);

    ptrdiff_t lb;
    ptrdiff_t send_elsize;
    ptrdiff_t recv_elsize;

    MPI_Type_get_extent(sendtype, &lb, &send_elsize);
    MPI_Type_get_extent(recvtype, &lb, &recv_elsize);

#ifndef NO_ISEND_IRECV_IN_DOMAIN
    int n_requests;
    MPI_Request requests[NTask * 2];
    n_requests = 0;


    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;

        if(target >= NTask) continue;
        if(recvcnts[target] == 0) continue;
        MPI_Irecv(
                ((char*) recvbuf) + recv_elsize * rdispls[target], 
                recvcnts[target],
                recvtype, target, 101934, comm, &requests[n_requests++]);
    }

    MPI_Barrier(comm);	/* not really necessary, but this will guarantee that all receives are
                                       posted before the sends, which helps the stability of MPI on 
                                       bluegene, and perhaps some mpich1-clusters */

    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;
        if(target >= NTask) continue;
        if(sendcnts[target] == 0) continue;
        MPI_Isend(((char*) sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target],
                sendtype, target, 101934, comm, &requests[n_requests++]);
    }

    MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

#else
    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;

        if(target >= NTask) continue;
        if(sendcnts[target] == 0 && recvcnts[target] == 0) continue;
        MPI_Sendrecv(((char*)sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target], sendtype, 
                target, 101934,
                ((char*)recvbuf) + recv_elsize * rdispls[target],
                recvcnts[target], recvtype, 
                target, 101934, 
                comm, MPI_STATUS_IGNORE);

    }
#endif
    /* ensure the collective-ness */
    MPI_Barrier(comm);

    return 0;
} 
