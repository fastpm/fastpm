#include <string.h>
#include <mpi.h>

#include <fftw3.h>
#include <fftw3-mpi.h>
#include <pfft.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/transfer.h>

#include "pmpfft.h"
static MPI_Datatype MPI_PTRDIFF = (MPI_Datatype) 0;

#if FASTPM_FFT_PRECISION == 64
    #define plan_dft_r2c pfft_plan_dft_r2c
    #define plan_dft_c2r pfft_plan_dft_c2r
    #define execute_dft_r2c pfft_execute_dft_r2c
    #define execute_dft_c2r pfft_execute_dft_c2r
    #define plan_dft_r2c_fftw fftw_mpi_plan_dft_r2c
    #define plan_dft_c2r_fftw fftw_mpi_plan_dft_c2r
    #define execute_dft_r2c_fftw fftw_mpi_execute_dft_r2c
    #define execute_dft_c2r_fftw fftw_mpi_execute_dft_c2r
    #define _pfft_init pfft_init
    #define _pfft_cleanup pfft_cleanup
    #define destroy_plan pfft_destroy_plan
    #define destroy_plan_fftw fftw_destroy_plan

#elif FASTPM_FFT_PRECISION == 32
    #define plan_dft_r2c pfftf_plan_dft_r2c
    #define plan_dft_c2r pfftf_plan_dft_c2r
    #define plan_dft_r2c_fftw fftwf_mpi_plan_dft_r2c
    #define plan_dft_c2r_fftw fftwf_mpi_plan_dft_c2r
    #define execute_dft_r2c pfftf_execute_dft_r2c
    #define execute_dft_c2r pfftf_execute_dft_c2r
    #define execute_dft_r2c_fftw fftwf_mpi_execute_dft_r2c
    #define execute_dft_c2r_fftw fftwf_mpi_execute_dft_c2r
    #define _pfft_init pfftf_init
    #define _pfft_cleanup pfftf_cleanup
    #define destroy_plan pfftf_destroy_plan
    #define destroy_plan_fftw fftwf_destroy_plan
#endif

void
pm_module_init() 
{
    if(MPI_PTRDIFF) return;
        
    _pfft_init();

    if(sizeof(ptrdiff_t) == 8) {
        MPI_PTRDIFF = MPI_LONG;
    } else {
        MPI_PTRDIFF = MPI_INT;
    }
}

void 
pm_module_cleanup() 
{
    if(!MPI_PTRDIFF) return;
    _pfft_cleanup();

    MPI_PTRDIFF = (MPI_Datatype) 0;
}

static size_t fftw_local_size_dft_r2c(int nrnk, ptrdiff_t * n, MPI_Comm comm,
                        int flags, 
                        ptrdiff_t * isize, ptrdiff_t * istart,
                        ptrdiff_t * osize, ptrdiff_t * ostart
                        ) {
    size_t allocsize;
    ptrdiff_t n2[nrnk];
    int i;
    for(i = 0; i < nrnk; i++ ) {
        n2[i] = n[i];
    }
    /* r2c is always padded !*/
    n2[nrnk - 1] = n[nrnk - 1] / 2 + 1;

    /* translate to a compatible interface with PFFT */
    for(i = 0; i < nrnk; i ++ ){
        istart[i] = 0;
        ostart[i] = 0;
        isize[i] = n2[i];
        if(i == nrnk - 1) {
            /* real input */
            isize[i] *= 2;
        }
        osize[i] = n2[i];
    }

    if(FFTW_MPI_TRANSPOSED_OUT & flags) {
        allocsize = fftw_mpi_local_size_transposed(
                3, n2, comm, 
                &isize[0], &istart[0],
                &osize[1], &ostart[1]);
    } else {
        allocsize = fftw_mpi_local_size(
                3, n2, comm, 
                &isize[0], &istart[0]);
        osize[0] = isize[0];
        osize[1] = isize[1];
    }
    return allocsize;
}

void pm_init(PM * pm, PMInit * init, MPI_Comm comm) {

    pm->init = *init;
    pm->mem = _libfastpm_get_gmem();

    /* initialize the domain */
    MPI_Comm_rank(comm, &pm->ThisTask);
    MPI_Comm_size(comm, &pm->NTask);

    int Ny = init->NprocY;
    int Nx;
    if(Ny <= 0) {
        Ny = 1;
        Nx = pm->NTask;
        if(!init->use_fftw) {
            for(; Ny * Ny < pm->NTask; Ny ++) continue;
            for(; Ny >= 1; Ny--) {
                if (pm->NTask % Ny == 0) break;
                continue;
            }
        }
    } else {
        if(pm->NTask % Ny != 0) {
            fastpm_raise(-1, "NprocY(%d) and NTask(%d) is incompatible\n", Ny, pm->NTask);
        }
    }
    Nx = pm->NTask / Ny;
    pm->Nproc[0] = Nx;
    pm->Nproc[1] = Ny;
    if(init->use_fftw) {
        if(Ny != 1) {
            fastpm_raise(-1, "FFTW requires Ny == 1; Ny = %d\n", Ny);
        }
    }
    int d;
    if(init->Nmesh % 2 != 0) {
        fastpm_raise(-1, "Nmesh must be even, but %d is odd.\n", init->Nmesh);
    }
    pm->Norm = 1.0;
    pm->Volume = 1.0;
    for(d = 0; d < 3; d ++) {
        pm->Nmesh[d] = init->Nmesh;
        pm->BoxSize[d] = init->BoxSize;

        pm->CellSize[d] = pm->BoxSize[d] / pm->Nmesh[d];
        pm->InvCellSize[d] = 1.0 / pm->CellSize[d]; 
        pm->Norm *= pm->Nmesh[d];
        pm->Volume *= pm->BoxSize[d];

    }


    pfft_create_procmesh(2, comm, pm->Nproc, &pm->Comm2D);

    if(init->use_fftw) {
        pm->allocsize = 2 * fftw_local_size_dft_r2c(
                3, pm->Nmesh, pm->Comm2D, 
                (pm->init.transposed?FFTW_MPI_TRANSPOSED_OUT:0),
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);
    } else {
        pm->allocsize = 2 * pfft_local_size_dft_r2c(
                3, pm->Nmesh, pm->Comm2D, 
                (pm->init.transposed?PFFT_TRANSPOSED_OUT:0)
                | PFFT_PADDED_R2C, 
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);
    }

    /* Note that we need to fix up the padded size of the real data;
     * and transpose with strides , */


    pm->IRegion.strides[2] = 1;
    pm->IRegion.strides[1] = pm->IRegion.size[2];
    pm->IRegion.strides[0] = pm->IRegion.size[1] * pm->IRegion.strides[1];
    pm->IRegion.total = pm->IRegion.size[0] * pm->IRegion.strides[0];

    /* remove padding from the view */
    pm->IRegion.size[2] = pm->Nmesh[2];

    if(pm->init.transposed) {
        if(pm->init.use_fftw) {
            /* FFTW transposed, y, x, z */
            pm->ORegion.strides[2] = 1;
            pm->ORegion.strides[0] = pm->ORegion.size[2];
            pm->ORegion.strides[1] = pm->ORegion.size[0] * pm->ORegion.strides[0];
            pm->ORegion.total = pm->ORegion.size[1] * pm->ORegion.strides[1];

        } else {
            /* PFFT transposed, y, z, x */
            pm->ORegion.strides[0] = 1;
            pm->ORegion.strides[2] = pm->ORegion.size[0];
            pm->ORegion.strides[1] = pm->ORegion.size[2] * pm->ORegion.strides[2];
            pm->ORegion.total = pm->ORegion.size[1] * pm->ORegion.strides[1];
        }
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

        /* Here we sum the size instead of using the offset,
         * because pfft reports local_i_start and local_ni of zero.
         * if the rank does not contain any blocks. Using the offsets
         * will cause non-increasing edges.*/

        MPI_Cart_sub(pm->Comm2D, remain_dims, &projected);
        MPI_Allgather(&pm->IRegion.size[d], 1, MPI_PTRDIFF, 
            pm->Grid.edges_int[d], 1, MPI_PTRDIFF, projected);

        MPI_Comm_free(&projected);
        int j;
        ptrdiff_t sum = 0;
        for(j = 0; j <= pm->Nproc[d]; j ++) {
            ptrdiff_t sum1 = sum + pm->Grid.edges_int[d][j];
            pm->Grid.edges_int[d][j] = sum;
            sum = sum1;
        }

        for(j = 0; j <= pm->Nproc[d]; j ++) {
            pm->Grid.edges_float[d][j] = 1.0 * pm->Grid.edges_int[d][j] / pm->Nmesh[d] * pm->BoxSize[d];
        }

        /* Last edge is at the edge of the box */
        int last = pm->Nproc[d];
        pm->Grid.edges_float[d][last] = pm->BoxSize[d];
        if (pm->Grid.edges_int[d][last] != pm->Nmesh[d]) {
            fastpm_raise(-1, "last edge did not match Nmesh. Internal Error.\n");

        };
        /* fill in the look up table */
        for(j = 0; j < pm->Nproc[d]; j ++) {
            int i;
            for(i = pm->Grid.edges_int[d][j]; i < pm->Grid.edges_int[d][j+1]; i ++) {
                pm->Grid.MeshtoCart[d][i] = j;
            }
        }
    }

    FastPMFloat * canvas = pm_alloc(pm);
    FastPMFloat * workspace = pm_alloc(pm);

    if(pm->init.use_fftw) {
        pm->r2c = plan_dft_r2c_fftw(
                3, pm->Nmesh, (void*) workspace, (void*) canvas, 
                pm->Comm2D, 
                (pm->init.transposed?FFTW_MPI_TRANSPOSED_OUT:0)
                | FFTW_ESTIMATE 
                | FFTW_DESTROY_INPUT
                );
        pm->c2r = plan_dft_c2r_fftw(
                3, pm->Nmesh, (void*) canvas, (void*) canvas, 
                pm->Comm2D, 
                (pm->init.transposed?FFTW_MPI_TRANSPOSED_IN:0)
                | FFTW_ESTIMATE 
                | FFTW_DESTROY_INPUT
                );
    } else {
        pm->r2c = plan_dft_r2c(
                3, pm->Nmesh, (void*) workspace, (void*) canvas, 
                pm->Comm2D,
                PFFT_FORWARD, 
                (pm->init.transposed?PFFT_TRANSPOSED_OUT:0)
                | PFFT_PADDED_R2C 
                | PFFT_ESTIMATE 
                | PFFT_TUNE
                //| PFFT_MEASURE
                | PFFT_DESTROY_INPUT
                );
        pm->c2r = plan_dft_c2r(
                3, pm->Nmesh, (void*) workspace, (void*) workspace, 
                pm->Comm2D,
                PFFT_BACKWARD, 
                (pm->init.transposed?PFFT_TRANSPOSED_IN:0)
                | PFFT_PADDED_C2R 
                | PFFT_ESTIMATE 
                //| PFFT_MEASURE
                | PFFT_TUNE
                | PFFT_DESTROY_INPUT
                );
    }

    pm_free(pm, workspace);
    pm_free(pm, canvas);

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

void 
pm_destroy(PM * pm) 
{
    int d;
    if(pm->init.use_fftw) {
        destroy_plan_fftw(pm->r2c);
        destroy_plan_fftw(pm->c2r);
    } else {
        destroy_plan(pm->r2c);
        destroy_plan(pm->c2r);
    }
    for(d = 0; d < 3; d++) {
        free(pm->MeshtoK[d]);
    }

    for(d = 0; d < 2; d++) {
        free(pm->Grid.MeshtoCart[d]);
        free(pm->Grid.edges_int[d]);
        free(pm->Grid.edges_float[d]);
    }
    MPI_Comm_free(&pm->Comm2D);
}

int pm_pos_to_rank(PM * pm, double pos[3]) {
    int d;
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
        if(UNLIKELY(ipos < 0)) {
            ipos = ipos % pm->Nmesh[d];
            if(ipos < 0) { ipos += pm->Nmesh[0]; }
        }
        if(UNLIKELY(ipos >= pm->Nmesh[d])) {
            ipos = ipos % pm->Nmesh[d];
        }
        rank2d[d] = pm->Grid.MeshtoCart[d][ipos];
    }
    return rank2d[0] * pm->Nproc[1] + rank2d[1];
}

void pm_r2c(PM * pm, FastPMFloat * from, FastPMFloat * to) {
    VALGRIND_CHECK_MEM_IS_DEFINED(from, sizeof(from[0]) * pm->allocsize);

    /* A gaussian of variance 1 becomes a complex gausian of variance 1/2 * (1 / Norm) in real and imag */

    /* workspace to canvas*/
    if(pm->init.use_fftw) {
        execute_dft_r2c_fftw(pm->r2c, from, (void*)to);
    } else {
        execute_dft_r2c(pm->r2c, from, (void*)to);
    }
    ptrdiff_t i;
#pragma omp parallel for
    for(i = 0; i < pm->allocsize; i ++) {
        to[i] *= 1 / pm->Norm;
    }

    VALGRIND_MAKE_MEM_DEFINED(to, sizeof(to[0]) * pm->allocsize);
}

void pm_c2r(PM * pm, FastPMFloat * inplace) {
    /* r2c and c2r round trip is unitary */
    VALGRIND_CHECK_MEM_IS_DEFINED(inplace, sizeof(inplace[0]) * pm->allocsize);
    if(pm->init.use_fftw) {
        execute_dft_c2r_fftw(pm->c2r, (void*) inplace, inplace);
    } else {
        execute_dft_c2r(pm->c2r, (void*) inplace, inplace);
    }
    VALGRIND_MAKE_MEM_DEFINED(inplace, sizeof(inplace[0]) * pm->allocsize);
}

#define unravel(ind, i, d0, d1, d2, strides) \
i[d0] = ind / strides[d0]; ind %= strides[d0]; \
i[d1] = ind / strides[d1]; ind %= strides[d1]; \
i[d2] = ind

void pm_unravel_o_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]) {
    /*
     * using pm_unravel_o_index function is slower than pm_inc_o_index, thus it is only used
     * during dev to test pm_inc_o_index.
     * */
    ptrdiff_t tmp = ind;
    if(pm->init.transposed) {
        if(pm->init.use_fftw) {
            /* y, x, z*/
            unravel(tmp, i, 1, 0, 2, pm->ORegion.strides);
        } else {
            /* y, z, x*/
            unravel(tmp, i, 1, 2, 0, pm->ORegion.strides);
        }
    } else {
        unravel(tmp, i, 0, 1, 2, pm->ORegion.strides);
    }
}
void pm_unravel_i_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]) {
    ptrdiff_t tmp = ind;
    unravel(tmp, i, 0, 1, 2, pm->IRegion.strides);
}

ptrdiff_t pm_ravel_o_index(PM * pm, ptrdiff_t i[3]) {
    ptrdiff_t ind = 0;
    int d;
    for(d = 0; d < 3; d++) {
        ind += pm->ORegion.strides[d] * i[d];
    }
    return ind;
}
ptrdiff_t pm_ravel_i_index(PM * pm, ptrdiff_t i[3]) {
    ptrdiff_t ind = 0;
    int d;
    for(d = 0; d < 3; d++) {
        ind += pm->IRegion.strides[d] * i[d];
    }
    return ind;
}

#define inc(i, d0, d1, d2, size) \
            i[d2] ++; \
            if(UNLIKELY(i[d2] == size[d2])) { \
                i[d2] = 0; i[d1] ++; \
                if(UNLIKELY(i[d1] == size[d1])) { \
                    i[d1] = 0; \
                    i[d0] ++; \
                } \
            }

/* returns number of items moved in linear index */
int pm_inc_o_index(PM * pm, ptrdiff_t i[3]) {
    if(pm->init.transposed) {
        if(pm->init.use_fftw) {
            /* y, x, z */
            inc(i, 1, 0, 2, pm->ORegion.size);
        } else {
            /* y, z, x */
            inc(i, 1, 2, 0, pm->ORegion.size);
        }
    } else {
        /* x, y, z*/
        inc(i, 0, 1, 2, pm->ORegion.size);
    }
    return 1;
}

/* returns number of items moved in linear index */
int pm_inc_i_index(PM * pm, ptrdiff_t i[3]) {
    /* can't use the macro because of the padding */
    int rt = 1;
    i[2] ++;
    if(UNLIKELY(i[2] == pm->IRegion.size[2])) { /* the padding !*/
        rt += 2;
        i[2] = 0;
        i[1] ++;
        if(UNLIKELY(i[1] == pm->IRegion.size[1])) {
            i[1] = 0;
            i[0] ++;
        }
    }
    return rt;
}

int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {

    int ThisTask;
    int NTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    {
        /* if the send is dense, use MPI_Alltoallv directly. */
        int i;

        size_t send_requests = 0;

        for(i = 0; i < NTask; i ++) {
            if(sendcnts[i] > 0) {
                send_requests ++;
            }
        }
        int dense = send_requests > 128;
        MPI_Allreduce(MPI_IN_PLACE, &dense, 1, MPI_INT, MPI_SUM, comm);
        /* dense is number of ranks does a lot of sends. */
        if (dense > 0) {
            fastpm_info("Using MPI's Alltoallv");
            return MPI_Alltoallv(
                sendbuf, sendcnts, sdispls, sendtype,
                recvbuf, recvcnts, rdispls, recvtype,
                comm);
        } /* else */
    }
    fastpm_info("Using sparse Alltoallv");

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

        VALGRIND_CHECK_MEM_IS_DEFINED(sendbuf, send_elsize * sdispls[target]);

        MPI_Isend(((char*) sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target],
                sendtype, target, 101934, comm, &requests[n_requests++]);
    }

    MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

    for(ngrp = 0; ngrp < (1 << PTask); ngrp++) {
        int target = ThisTask ^ ngrp;
        if(target >= NTask) continue;
        if(recvcnts[target] == 0) continue;

        VALGRIND_MAKE_MEM_DEFINED(recvbuf, recv_elsize * rdispls[target]);
    }
#else
    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;

        if(target >= NTask) continue;
        if(sendcnts[target] == 0 && recvcnts[target] == 0) continue;

        VALGRIND_CHECK_MEM_IS_DEFINED(sendbuf, send_elsize * sdispls[target]);

        MPI_Sendrecv(((char*)sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target], sendtype, 
                target, 101934,
                ((char*)recvbuf) + recv_elsize * rdispls[target],
                recvcnts[target], recvtype, 
                target, 101934, 
                comm, MPI_STATUS_IGNORE);

        VALGRIND_MAKE_MEM_DEFINED(recvbuf, recv_elsize * rdispls[target]);
    }
#endif
    /* ensure the collective-ness */
    MPI_Barrier(comm);

    return 0;
} 

