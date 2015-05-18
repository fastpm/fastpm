/**
 *  2D Domain Decomposition
 *
 *  This code builds 2d pencel domains for PM simulations
 *
 *  Currently we use FFTW thus the domains are reduced to slabs.
 * 
 *  Particles near the edges are ghosts.
 *
 *  Author: Yu Feng <rainwoodman@gmail.com>
 */

#define MAKE_ASSERTIONS

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <assert.h>
#include <math.h>
#include "msg.h"
#include "particle.h"
#include <fftw3-mpi.h>

#include "heap.h"
#include "msort.c"

double BoxSize;
int * MeshToTask[2];
int NTask;
int ThisTask;
int NTask2D[2];
int TaskStrides[2];
ptrdiff_t local_n[2];
ptrdiff_t local_start[2];

MPI_Comm Comm2D;
ptrdiff_t Nmesh;
double MeshPerBox;

void domain_init(int Nmesh_, double BoxSize_) {
    Nmesh = Nmesh_;
    BoxSize = BoxSize_;
    MeshPerBox = Nmesh/ BoxSize;

    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    /* fake it with slabs, will be pfft at some point */

    NTask2D[0] = NTask;
    NTask2D[1] = 1;
    TaskStrides[0] = 1; 
    TaskStrides[1] = 1; 
    MPI_Cart_create(MPI_COMM_WORLD, 2, NTask2D, NTask2D, 0, &Comm2D);

    fftwf_mpi_local_size_3d(Nmesh, Nmesh, Nmesh, MPI_COMM_WORLD,
			  &local_n[0], &local_start[0]);
    local_n[1] = Nmesh;
    local_start[1] = 0;
    
    /* build array for X and Y */
    for(int d = 0; d < 2; d ++) {
        MeshToTask[d] = malloc(sizeof(int) * Nmesh);
        MPI_Comm RANKS;
        int keep[2] = {0};
        keep[d] = 1; 
        MPI_Cart_sub(Comm2D, keep, &RANKS);

        int myrank;
        MPI_Comm_rank(RANKS, &myrank);
        int * recvcounts = alloca(sizeof(int) * NTask2D[d]);
        int * recvdispls = alloca(sizeof(int) * NTask2D[d]);
        
        int itmp = local_n[d];
        MPI_Allgather (&itmp, 1, MPI_INT, recvcounts, 1, MPI_INT, RANKS);
        int * buf = alloca(sizeof(int) * local_n[d]);
        for(int i = 0; i < local_n[d]; i ++) {
            buf[i] = myrank;
        }
        recvdispls[0] = 0;
        for(int i = 1; i < NTask2D[d]; i ++){
            recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
        }
        MPI_Allgatherv (buf, local_n[d], MPI_INT, 
                MeshToTask[d], recvcounts, recvdispls, MPI_INT, RANKS);
        MPI_Comm_free(&RANKS);
    }

}
void domain_finalize() {
    MPI_Comm_free(&Comm2D);
    free(MeshToTask[0]);
    free(MeshToTask[1]);
}
static int par_node(float x[3]) {
    int task = 0;
    for(int d = 0; d < 2; d ++) {
        int ind = floor(x[d] * MeshPerBox);
        while(ind < 0) ind += Nmesh;
        while(ind >= Nmesh) ind -= Nmesh;
        task += TaskStrides[d] * MeshToTask[d][ind];
    }
    return task;
}

/* this will sort my particles to the top, 
 * then the other particles near the end by order */
static int cmp_par_task(void * c1, void * c2, void * data) {
    int64_t i1 = *((int64_t *) c1);
    int64_t i2 = *((int64_t *) c2);
    Particles * p = data;
    int task1 = par_node(p->x[i1]);
    int task2 = par_node(p->x[i2]);
    if (task1 == ThisTask) task1 = -1;
    if (task2 == ThisTask) task2 = -1;
    int rt = (task1 > task2) - (task2 > task1);
    if(rt) return rt;
    long long id1 = p->id[i1];
    long long id2 = p->id[i2];
    return (id1 > id2) - (id2 > id1);
}
struct GhostBuf {
    float x[3];
    union {
        int TargetTask;
        int OriginalTask;
    };
    int index;
};

static int cmp_ghostbuf_target(void * c1, void * c2, void* data) {
    int t1 = ((struct GhostBuf*) c1)->TargetTask;
    int t2 = ((struct GhostBuf*) c2)->TargetTask;
    return (t1 > t2) - (t2 > t1);
}

int domain_create_ghosts(Particles* p, double eps) {

    int * SendCount = alloca(sizeof(int) * NTask);
    int * RecvCount = alloca(sizeof(int) * NTask);
    int * SendDispl = alloca(sizeof(int) * NTask);
    int * RecvDispl = alloca(sizeof(int) * NTask);

    /* FIXME: for now at most create np_local ghosts; we are CHECK this! */
    struct GhostBuf * sendbuf = heap_allocate(sizeof(struct GhostBuf) * p->np_local);

    for(int i = 0; i < NTask; i ++) {
        SendCount[i] = 0;
        RecvCount[i] = 0;
        SendDispl[i] = 0;
        RecvDispl[i] = 0;
    }

    /* build the send buf */
    int nsend = 0;
    for(int i = 0; i < p->np_local; i ++) {
        float * x = p->x[i];
        float tmp[3];
        memcpy(tmp, x, 3 * sizeof(float));
        /* check all neighbours to see if the particle
         * crosses the edge boundary */
        int exportedTasks[4] = {-1, -1, -1, -1};
        int exported = 0; 
        for(int xfac = 0; xfac <=1; xfac ++) 
        for(int yfac = 0; yfac <=1; yfac ++) {
            tmp[0] = x[0] + xfac * eps;
            tmp[1] = x[1] + yfac * eps;
            int task = par_node(tmp);
            if(task == ThisTask) goto skip;
            for(int j = 0; j < exported; j ++) {
                if(exportedTasks[j] == task) goto skip;
            }
            exportedTasks[exported++] = task;
            float * x1 = sendbuf[nsend].x;
            memcpy(x1, x, sizeof(float) * 3);
            sendbuf[nsend].TargetTask = task;
            sendbuf[nsend].index = i;
            SendCount[task] ++;
            nsend ++;
        skip:
            continue;
        }
    }
    size_t tmp_size = 2 * nsend * sizeof(void*) + sizeof(sendbuf[0]);
    void * tmp = heap_allocate(tmp_size);
    qsort_r_with_tmp(sendbuf, nsend, sizeof(sendbuf[0]), cmp_ghostbuf_target, p, tmp, tmp_size);
    heap_return(tmp);
    for(int i = 0; i < nsend; i ++) {
        sendbuf[i].OriginalTask = ThisTask;
    }

    /* now exchange the ghosts */
    struct GhostBuf * recvbuf = &sendbuf[nsend];

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    int nrecv = RecvDispl[NTask - 1] + RecvCount[NTask - 1];
    if( p->np_local + nrecv  > p->np_allocated) {
        msg_abort(2415, "Not enough particle space ghosts, nrecv=%d, left=%td",
            nrecv, p->np_allocated - p->np_local
            );
    }

    MPI_Datatype MPI_SENDBUF;
    MPI_Type_contiguous(sizeof(sendbuf[0]), MPI_BYTE, &MPI_SENDBUF);
    MPI_Type_commit(&MPI_SENDBUF);

    MPI_Alltoallv(sendbuf, SendCount, SendDispl, MPI_SENDBUF,
                  recvbuf, RecvCount, RecvDispl, MPI_SENDBUF, 
                  MPI_COMM_WORLD);

    MPI_Type_free(&MPI_SENDBUF);
    
    /* add ghosts to particles */
    for(int i = 0; i < nrecv; i ++) {
        memcpy(p->x[p->np_local + i], recvbuf[i].x, sizeof(float) * 3);
        /* ghost ID encodes the original index and task */
        p->id[p->np_local + i] = (((int64_t)recvbuf[i].OriginalTask) << 32L) + recvbuf[i].index;
    }
    int maxghosts = 0;
    MPI_Allreduce(&nrecv, &maxghosts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    msg_printf(verbose, "Max number of ghosts is %d, %g / 1\n", maxghosts, 1.0 * maxghosts/ p->np_allocated);

    heap_return(sendbuf);
    return nrecv;
}

void domain_annihilate_ghosts(Particles * p, 
        int nghosts,
        float (* f3)[3])
{
    int np_local = p->np_local;
    int * SendCount = alloca(sizeof(int) * NTask);
    int * RecvCount = alloca(sizeof(int) * NTask);
    int * SendDispl = alloca(sizeof(int) * NTask);
    int * RecvDispl = alloca(sizeof(int) * NTask);

    for(int i = 0; i < NTask; i ++) {
        SendCount[i] = 0;
        RecvCount[i] = 0;
        SendDispl[i] = 0;
        RecvDispl[i] = 0;
    }

    float (*f3g)[3] = f3 + np_local;

    int * indexsend = heap_allocate(sizeof(int*) * nghosts);

    for(int i = 0; i < nghosts; i ++) {
        int originalIndex = p->id[p->np_local + i] & 0xffffffff;
        int originalTask = p->id[p->np_local + i] >> 32L;
        SendCount[originalTask] ++;
        indexsend[i] = originalIndex;
    }

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    int nrecv = RecvDispl[NTask - 1] + RecvCount[NTask - 1];

    int * indexrecv = heap_allocate(sizeof(int*) * nrecv);

    MPI_Alltoallv(indexsend, SendCount, SendDispl, MPI_INT,
                indexrecv, RecvCount, RecvDispl, MPI_INT,
                MPI_COMM_WORLD);

    /*scratch space */
    float (*g3)[3] = heap_allocate(sizeof(float) * 3 * nrecv);

    MPI_Datatype MPI_F3;
    MPI_Type_contiguous(3, MPI_FLOAT, &MPI_F3);
    MPI_Type_commit(&MPI_F3);
    MPI_Alltoallv(f3g, SendCount, SendDispl, MPI_F3,
                g3, RecvCount, RecvDispl, MPI_F3,
                MPI_COMM_WORLD);
    MPI_Type_free(&MPI_F3);

#ifdef MAKE_ASSERTIONS
    double sum0[3] = {0};
    for(int i = 0; i < np_local + nghosts; i ++) {
        for(int d = 0; d < 3; d ++) {
            sum0[d] += f3[i][d];
        }
    }
#endif

    for(int i = 0; i < nrecv; i ++) {
        int index = indexrecv[i];
        if(index >= np_local) {
            msg_abort(9934, "Bad index");
        }
        for(int d = 0; d < 3; d ++) {
            f3[index][d] += g3[i][d];
        }
    }
#ifdef MAKE_ASSERTIONS
    double sum1[3] = {0};
    for(int i = 0; i < np_local; i ++) {
        for(int d = 0; d < 3; d ++) {
            sum1[d] += f3[i][d];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Checking ghost forces: %g %g %g ~ %g %g %g \n", 
            sum0[0], sum0[1], sum0[2], 
            sum1[0], sum0[1], sum0[2]);
#endif
    heap_return(g3);
    heap_return(indexrecv);
    heap_return(indexsend);
}

void domain_wrap(Particles * p) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < p->np_local; i ++) {
        for(int d = 0; d < 3; d++) {
            while(p->x[i][d] < 0) p->x[i][d] += BoxSize;
            while(p->x[i][d] >= BoxSize) p->x[i][d] -= BoxSize;
        }
    }
}
static void permute(void * data, int np, size_t elsize, int64_t * ind) {
    void * tmp = heap_allocate(np * elsize);
    char * p = tmp;
    char * q = data;
    for(int i = 0; i < np; i ++) {
        memcpy(&p[i * elsize], &q[ind[i] * elsize], elsize);
    }
    memcpy(data, tmp, elsize * np);
    heap_return(tmp);
}
struct DomainBuf {
    float x[3];
    float v[3];
    float dx1[3];
    float dx2[3];
    int64_t id;
};

void domain_decompose(Particles* p) {
#ifdef MAKE_ASSERTIONS
    {
        size_t total1 = p->np_local;
        size_t total = 0;
        MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        msg_printf(verbose, "Total number of particles = %td\n", total);
        total1 = 0;
        for(int i = 0; i < p->np_local; i ++) {
            total1 += p->id[i];
        }
        MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        msg_printf(verbose, "Sum of particle ID = %td\n", total);
    }
#endif
    int64_t * ind = heap_allocate(sizeof(int64_t) * p->np_local);
    for(int i = 0; i < p->np_local; i++) {
        ind[i] = i;
    }
    size_t tmp_size = 2 * p->np_local * sizeof(void*) + sizeof(int64_t);
    void * tmp = heap_allocate(tmp_size);
    qsort_r_with_tmp(ind, p->np_local, sizeof(int64_t), cmp_par_task, p, tmp, tmp_size);
    heap_return(tmp);

    permute(p->x, p->np_local, sizeof(float) * 3, ind);
    permute(p->v, p->np_local, sizeof(float) * 3, ind);
    permute(p->id, p->np_local, sizeof(int64_t), ind);

    if(p->dx1)
        permute(p->dx1, p->np_local, sizeof(float) * 3, ind);
    if(p->dx2)
        permute(p->dx2, p->np_local, sizeof(float) * 3, ind);

    heap_return(ind);

    int * SendCount = alloca(sizeof(int) * NTask);
    int * RecvCount = alloca(sizeof(int) * NTask);
    int * SendDispl = alloca(sizeof(int) * NTask);
    int * RecvDispl = alloca(sizeof(int) * NTask);
    for(int i = 0; i < NTask; i ++) {
        SendCount[i] = 0;
        RecvCount[i] = 0;
        SendDispl[i] = 0;
        RecvDispl[i] = 0;
    }
    for(int i = 0; i < p->np_local; i ++) {
        int task = par_node(p->x[i]);
        SendCount[task] ++;
    }
    int nsend = p->np_local - SendCount[ThisTask];
    /* chop the p buffer to local and send buf */
    p->np_local = SendCount[ThisTask];
    /* local particles are not sent */
    SendCount[ThisTask] = 0;

    struct DomainBuf * sendbuf = heap_allocate(sizeof(struct DomainBuf) * nsend);

    for(int i = 0; i < nsend; i ++) {
        for(int d = 0; d < 3; d ++) {
            sendbuf[i].x[d] = p->x[p->np_local + i][d];
            sendbuf[i].v[d] = p->v[p->np_local + i][d];
            if(p->dx1)
                sendbuf[i].dx1[d] = p->dx1[p->np_local + i][d];
            if(p->dx2)
                sendbuf[i].dx2[d] = p->dx2[p->np_local + i][d];
        }
        sendbuf[i].id = p->id[p->np_local + i];
    }
    MPI_Datatype MPI_PAR;
    MPI_Type_contiguous(sizeof(struct DomainBuf), MPI_BYTE, &MPI_PAR);
    MPI_Type_commit(&MPI_PAR);

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    /* prepare the recv buf from scratch */
    size_t nrecv = RecvDispl[NTask - 1] + RecvCount[NTask - 1];

    struct DomainBuf * recvbuf = heap_allocate(sizeof(struct DomainBuf) * nrecv);

    MPI_Alltoallv(sendbuf, SendCount, SendDispl, MPI_PAR, 
            recvbuf, RecvCount, RecvDispl, MPI_PAR, MPI_COMM_WORLD);
    
    if(nrecv + p->np_local > p->np_allocated) {
        msg_abort(2314, "Not enough memory in particles for this exchange."
                "Need %d particles. Allocated %d particles\n",
                nrecv + p->np_local, p->np_allocated);

    }

    for(int i = 0; i < nrecv; i ++) {
        for(int d = 0; d < 3; d ++) {
            p->x[p->np_local][d] = recvbuf[i].x[d];
            p->v[p->np_local][d] = recvbuf[i].v[d];
            if(p->dx1)
                p->dx1[p->np_local][d] = recvbuf[i].dx1[d];
            if(p->dx2)
                p->dx2[p->np_local][d] = recvbuf[i].dx2[d];
        }
        p->id[p->np_local] = recvbuf[i].id;
        p->np_local ++;
    }

    heap_return(recvbuf);
    heap_return(sendbuf);

    MPI_Type_free(&MPI_PAR);

    for(int i = 0; i < p->np_local; i ++) {
        if(par_node(p->x[i]) != ThisTask) {
            msg_abort(2314, "Some particles are not in the correct domain after exchange.");
        }
    }
#ifdef MAKE_ASSERTIONS
    {
        size_t total1 = p->np_local;
        size_t total = 0;
        MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        msg_printf(verbose, "Total number of particles = %td\n", total);
        total1 = 0;
        for(int i = 0; i < p->np_local; i ++) {
            total1 += p->id[i];
        }
        MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        msg_printf(verbose, "Sum of particle ID = %td\n", total);
    }
#endif
}
