#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <omp.h>

#include "pmpfft.h"
#include "msg.h"
#include <fftw3.h>
#include <fftw3-mpi.h>

static MPI_Datatype MPI_PTRDIFF = NULL;

#if FFT_PRECISION == 64
    #define plan_dft_r2c pfft_plan_dft_r2c
    #define plan_dft_c2r pfft_plan_dft_c2r
    #define execute_dft_r2c pfft_execute_dft_r2c
    #define execute_dft_c2r pfft_execute_dft_c2r
    #define plan_dft_r2c_fftw fftw_mpi_plan_dft_r2c
    #define plan_dft_c2r_fftw fftw_mpi_plan_dft_c2r
    #define execute_dft_r2c_fftw fftw_mpi_execute_dft_r2c
    #define execute_dft_c2r_fftw fftw_mpi_execute_dft_c2r
    #define _pfft_init pfft_init
    #define destroy_plan pfft_destroy_plan
    #define destroy_plan_fftw fftw_destroy_plan

#elif FFT_PRECISION == 32
    #define plan_dft_r2c pfftf_plan_dft_r2c
    #define plan_dft_c2r pfftf_plan_dft_c2r
    #define plan_dft_r2c_fftw fftwf_mpi_plan_dft_r2c
    #define plan_dft_c2r_fftw fftwf_mpi_plan_dft_c2r
    #define execute_dft_r2c pfftf_execute_dft_r2c
    #define execute_dft_c2r pfftf_execute_dft_c2r
    #define execute_dft_r2c_fftw fftwf_mpi_execute_dft_r2c
    #define execute_dft_c2r_fftw fftwf_mpi_execute_dft_c2r
    #define _pfft_init pfftf_init
    #define destroy_plan pfftf_destroy_plan
    #define destroy_plan_fftw fftwf_destroy_plan
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

void pm_init(PM * pm, PMInit * init, PMIFace * iface, MPI_Comm comm) {

    module_init();

    pm->init = *init;
    pm->iface = *iface;

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
            msg_abort(-1, "NprocY(%d) and NTask(%d) is incompatible\n", Ny, pm->NTask);
        }
    }
    Nx = pm->NTask / Ny;
    pm->Nproc[0] = Nx;
    pm->Nproc[1] = Ny;
    if(init->use_fftw) {
        if(Ny != 1) {
            msg_abort(-1, "FFTW requires Ny == 1; Ny = %d\n", Ny);
        }
    }
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
    msg_printf(info, "ProcMesh : %d x %d ( %d Threads)\n", pm->Nproc[0], pm->Nproc[1], omp_get_max_threads());
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

    if(pm->init.use_fftw) {
        pm->r2c = plan_dft_r2c_fftw(
                3, pm->Nmesh, (void*) pm->workspace, (void*) pm->canvas, 
                pm->Comm2D, 
                (pm->init.transposed?FFTW_MPI_TRANSPOSED_OUT:0)
                | FFTW_ESTIMATE 
                | FFTW_DESTROY_INPUT
                );
        pm->c2r = plan_dft_c2r_fftw(
                3, pm->Nmesh, (void*) pm->workspace, (void*) pm->canvas, 
                pm->Comm2D, 
                (pm->init.transposed?FFTW_MPI_TRANSPOSED_IN:0)
                | FFTW_ESTIMATE 
                | FFTW_DESTROY_INPUT
                );
    } else {
        pm->r2c = plan_dft_r2c(
                3, pm->Nmesh, (void*) pm->workspace, (void*) pm->canvas, 
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
                3, pm->Nmesh, (void*) pm->workspace, (void*) pm->workspace, 
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

    pm->iface.free(pm->workspace);
    pm->iface.free(pm->canvas);
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

void 
pm_init_simple(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm) 
{
    PMInit pminit = {
        .Nmesh = Ngrid,
        .BoxSize = BoxSize,
        .NprocY = 0, /* 0 for auto, 1 for slabs */
        .transposed = 0,
        .use_fftw = 0,
    };

    pm_init(pm, &pminit, &p->iface, comm);

}

void 
pm_destroy(PM * pm) 
{
    pm_stop(pm);
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

void 
pm_start(PM * pm) 
{
    pm->canvas = pm->iface.malloc(sizeof(pm->canvas[0]) * pm->allocsize);
    pm->workspace = pm->iface.malloc(sizeof(pm->canvas[0]) * pm->allocsize);
}

void 
pm_stop(PM * pm) 
{
    if(pm->workspace)
        pm->iface.free(pm->workspace);
    if(pm->canvas && pm->canvas != pm->workspace)
        pm->iface.free(pm->canvas);

    pm->canvas = NULL;
    pm->workspace = NULL;
}

void pm_r2c(PM * pm) {
    if(pm->init.use_fftw) {
        execute_dft_r2c_fftw(pm->r2c, pm->workspace, (void*)pm->canvas);
    } else {
        execute_dft_r2c(pm->r2c, pm->workspace, (void*)pm->canvas);
    }
}

void pm_c2r(PM * pm) {
    if(pm->init.use_fftw) {
        execute_dft_c2r_fftw(pm->c2r, (void*) pm->workspace, pm->workspace);
    } else {
        execute_dft_c2r(pm->c2r, (void*) pm->workspace, pm->workspace);
    }
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

#define inc(i, d0, d1, d2, size) \
            i[d2] ++; \
            if(UNLIKELY(i[d2] == size[d2])) { \
                i[d2] = 0; i[d1] ++; \
                if(UNLIKELY(i[d1] == size[d1])) { \
                    i[d1] = 0; \
                    i[d0] ++; \
                } \
            }

void pm_inc_o_index(PM * pm, ptrdiff_t i[3]) {
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
}

void pm_inc_i_index(PM * pm, ptrdiff_t i[3]) {
    /* can't use the macro because of the padding */
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

static double 
sinc_unnormed(double x);

static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}

static double 
diff_kernel(double w) 
{
    /* order N = 1 super lanzcos kernel */
    /* 
     * This is the same as GADGET-2 but in fourier space: 
     * see gadget-2 paper and Hamming's book.
     * c1 = 2 / 3, c2 = 1 / 12
     * */
    return 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}

void 
pm_create_k_factors(PM * pm, PMKFactors * fac[3]) 
{ 
    /* This function populates fac with precalculated values that
     * are useful for force calculation. 
     * e.g. k**2 and the finite differentiation kernels. 
     * precalculating them means in the true kernel we only need a 
     * table look up. watch out for the offset ORegion.start
     * */
    int d;
    ptrdiff_t ind;
    for(d = 0; d < 3; d++) {
        fac[d] = malloc(sizeof(fac[0][0]) * pm->Nmesh[d]);
        double CellSize = pm->BoxSize[d] / pm->Nmesh[d];
        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            float k = pm->MeshtoK[d][ind];
            float w = k * CellSize;
            float ff = sinc_unnormed(0.5 * w);

            fac[d][ind].k_finite = 1 / CellSize * diff_kernel(w);
            fac[d][ind].kk_finite = k * k * ff * ff;
            fac[d][ind].kk = k * k;
            double tmp = sin(0.5 * k * CellSize);
            fac[d][ind].cic = 1 - 2. / 3 * tmp * tmp;
        }
    } 
}

void 
pm_destroy_k_factors(PM * pm, PMKFactors * fac[3]) 
{
    int d;
    for(d = 0; d < 3; d ++) {
        free(fac[d]);
    }
}

void 
pm_get_times(int istep,
    double time_step[],
    int nstep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1) 
{

    /* The last step is the terminal step. */
    *a_x = time_step[(istep >= nstep)?(nstep - 1):istep];
    *a_x1 = time_step[(istep + 1 >= nstep)?(nstep - 1):(istep + 1)];
    double a_xm1 = time_step[(istep > 0)?(istep - 1):0];
    *a_v = sqrt(a_xm1 * *(a_x));
    *a_v1 = sqrt(*a_x * *a_x1);
}

/* These shall go to another file. gravity stuff. */
#include "walltime.h"
static void 
apply_force_kernel(PM * pm, int dir) 
{
    /* This is the force in fourier space. - i k[dir] / k2 */

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k_finite;
            double kk_finite = 0;
            double kk = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
            }
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                pm->workspace[ind + 0] =   pm->canvas[ind + 1] * (k_finite / kk_finite);
                pm->workspace[ind + 1] = - pm->canvas[ind + 0] * (k_finite / kk_finite);
            } else {
                pm->workspace[ind + 0] = 0;
                pm->workspace[ind + 1] = 0;
            }
//            pm->workspace[ind + 0] = pm->canvas[ind + 0];
//            pm->workspace[ind + 1] = pm->canvas[ind + 1];
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);
}

void 
pm_calculate_forces(PMStore * p, PM * pm, double density_factor)
{
    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
        .nghosts = 0,
        .get_position = p->iface.get_position,
    };
    walltime_measure("/Force/Init");

    pm_append_ghosts(&pgd);
    walltime_measure("/Force/AppendGhosts");

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * we compensate this later at readout by density_factor.
     * */
    pm_paint(pm, p, p->np + pgd.nghosts);
    walltime_measure("/Force/Paint");
    
    pm_r2c(pm);
    walltime_measure("/Force/FFT");

#if 0
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("density-k.f4", "w"));
#endif

    /* calculate the forces save them to p->acc */

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        apply_force_kernel(pm, d);
        walltime_measure("/Force/Transfer");

#if 0
        char * fname[] = { "acc-0.f4", "acc-1.f4", "acc-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname[d], "w"));
#endif
        pm_c2r(pm);
        walltime_measure("/Force/FFT");

#if 0
{
    char buf[1000];
    sprintf(buf, "accr-%d.f4-rank-%d", d, pm->ThisTask);
    fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(buf, "w"));
}
#endif

#if 0
        char * fname2[] = { "accr-0.f4", "accr-1.f4", "accr-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname2[d], "w"));
#endif


#pragma omp parallel for
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            /* compensate the density is less than the true density */
            p->acc[i][d] = pm_readout_one(pm, p, i) * (density_factor / pm->Norm);
        }
        walltime_measure("/Force/Readout");

        pm_reduce_ghosts(&pgd, ACC[d]); 
        walltime_measure("/Force/ReduceGhosts");
    }
    pm_destroy_ghosts(&pgd);
    walltime_measure("/Force/Finish");

    MPI_Barrier(pm->Comm2D);
    walltime_measure("/Force/Wait");
}    

