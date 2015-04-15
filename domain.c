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
static int par_node(const float x[3]) {
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
static int cmp_par_task(const void * c1, const void * c2, void * data) {
    int task1 = par_node(((ParticleMinimum*) c1)->x);
    int task2 = par_node(((ParticleMinimum*) c2)->x);
    if (task1 == ThisTask) task1 = -1;
    if (task2 == ThisTask) task2 = -1;
    int rt = (task1 > task2) - (task2 > task1);
    if(rt) return rt;
    long long id1 = ((ParticleMinimum*) c1)->id;
    long long id2 = ((ParticleMinimum*) c2)->id;
    return (id1 > id2) - (id2 > id1);
}
struct SendBuf {
    float x[3];
    union {
        int TargetTask;
        int OriginalTask;
    };
    int index;
};
typedef struct {
    int * SendCount;    
    int * RecvCount;    
    int * SendDispl;    
    int * RecvDispl;    
} GhostData;

static int cmp_sendbuf_target(const void * c1, const void * c2, void* data) {
    int t1 = ((struct SendBuf*) c1)->TargetTask;
    int t2 = ((struct SendBuf*) c2)->TargetTask;
    return (t1 > t2) - (t2 > t1);
}

int domain_create_ghosts(Particles* const particles, double eps, void * scratch, size_t scratch_bytes)
{
    int * SendCount = alloca(sizeof(int) * NTask);
    int * RecvCount = alloca(sizeof(int) * NTask);
    int * SendDispl = alloca(sizeof(int) * NTask);
    int * RecvDispl = alloca(sizeof(int) * NTask);

    struct SendBuf * sendbuf = scratch;

    for(int i = 0; i < NTask; i ++) {
        SendCount[i] = 0;
        RecvCount[i] = 0;
        SendDispl[i] = 0;
        RecvDispl[i] = 0;
    }

    /* build the send buf */
    int nsend = 0;
    for(int i = 0; i < particles->np_local; i ++) {
        float * x = particles->p[i].x;
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
    void * tmp = &sendbuf[nsend];
    size_t tmp_size = scratch_bytes - sizeof(sendbuf[0]) * nsend;
    qsort_r_with_tmp(sendbuf, nsend, sizeof(sendbuf[0]), cmp_sendbuf_target, NULL, tmp, tmp_size);

    for(int i = 0; i < nsend; i ++) {
        sendbuf[i].OriginalTask = ThisTask;
    }

    /* now exchange the ghosts */
    struct SendBuf * recvbuf = &sendbuf[nsend];

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    int nrecv = RecvDispl[NTask - 1] + RecvCount[NTask - 1];
    if( particles->np_local + nrecv  > particles->np_allocated) {
        msg_abort(2415, "Not enough particle space ghosts, nrecv=%d, left=%td",
            nrecv, particles->np_allocated - particles->np_local
            );
    }
    if( (nrecv + nsend) * sizeof(sendbuf[0]) > scratch_bytes) {
        msg_abort(2415, "Not enough scratch space for building ghosts, nrecv=%d, nsend=%d, scratch=%d",
            nrecv, nsend, scratch_bytes / sizeof(sendbuf[0])
            );
    }

    MPI_Datatype MPI_SENDBUF;
    MPI_Type_contiguous(sizeof(sendbuf[0]), MPI_BYTE, &MPI_SENDBUF);
    MPI_Type_commit(&MPI_SENDBUF);

    MPI_Alltoallv(sendbuf, SendCount, SendDispl, MPI_SENDBUF,
                  recvbuf, RecvCount, RecvDispl, MPI_SENDBUF, 
                  MPI_COMM_WORLD);

    MPI_Type_free(&MPI_SENDBUF);
    
    Particle * ghosts = &particles->p[particles->np_local];
    /* add ghosts to particles */
    for(int i = 0; i < nrecv; i ++) {
        memset(&ghosts[i], 0, sizeof(ghosts[i]));
        memcpy(ghosts[i].x, recvbuf[i].x, sizeof(float) * 3);
        ghosts[i].OriginalTask = recvbuf[i].OriginalTask;
        ghosts[i].OriginalIndex = recvbuf[i].index;
        ghosts[i].id = -1;
    }
    int maxghosts = 0;
    MPI_Allreduce(&nrecv, &maxghosts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    msg_printf(verbose, "Max number of ghosts is %d, %g / 1\n", maxghosts, 1.0 * maxghosts/ particles->np_allocated);
    return nrecv;
}

void domain_annihilate_ghosts(Particles* const particles, 
        int nghosts,
        float (* f3)[3], 
        void * scratch, size_t scratch_bytes)
{
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
    char * base = scratch;

    /* count the number of ghosts to be returned */
    Particle * ghosts = &particles->p[particles->np_local];
    float (*f3g)[3] = f3 + particles->np_local;

    int * indexsend = (int*) base;
    base += sizeof(int) * nghosts;

    for(int i = 0; i < nghosts; i ++) {
        SendCount[ghosts[i].OriginalTask] ++;
        indexsend[i] = ghosts[i].OriginalIndex;
    }

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    int nrecv = RecvDispl[NTask - 1] + RecvCount[NTask - 1];

    int * indexrecv = (int*) base;
    base += sizeof(int) * nrecv;

    MPI_Alltoallv(indexsend, SendCount, SendDispl, MPI_INT,
                indexrecv, RecvCount, RecvDispl, MPI_INT,
                MPI_COMM_WORLD);

    /*scratch space */
    float (*g3)[3] = (float (*)[3]) base;

    MPI_Datatype MPI_F3;
    MPI_Type_contiguous(3, MPI_FLOAT, &MPI_F3);
    MPI_Type_commit(&MPI_F3);
    MPI_Alltoallv(f3g, SendCount, SendDispl, MPI_F3,
                g3, RecvCount, RecvDispl, MPI_F3,
                MPI_COMM_WORLD);
    MPI_Type_free(&MPI_F3);

    double sum0[3] = {0};
    for(int i = 0; i < particles->np_local + nghosts; i ++) {
        for(int d = 0; d < 3; d ++) {
            sum0[d] += f3[i][d];
        }
    }

    for(int i = 0; i < nrecv; i ++) {
        int index = indexrecv[i];
        if(index >= particles->np_local) {
            msg_abort(9934, "Bad index");
        }
        for(int d = 0; d < 3; d ++) {
            f3[index][d] += g3[i][d];
        }
    }
    double sum1[3] = {0};
    for(int i = 0; i < particles->np_local; i ++) {
        for(int d = 0; d < 3; d ++) {
            sum1[d] += f3[i][d];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Checking ghost forces: %g %g %g ~ %g %g %g \n", 
            sum0[0], sum0[1], sum0[2], 
            sum1[0], sum0[1], sum0[2]);
}

void domain_wrap(Particles * particles) {
    Particle * p = particles->p;
    for(int i = 0; i < particles->np_local; i ++) {
        for(int d = 0; d < 3; d++) {
            while(p[i].x[d] < 0) p[i].x[d] += BoxSize;
            while(p[i].x[d] >= BoxSize) p[i].x[d] -= BoxSize;
        }
    }
}
void domain_wrap_min(Snapshot * particles) {
    ParticleMinimum * p = particles->p;
    for(int i = 0; i < particles->np_local; i ++) {
        for(int d = 0; d < 3; d++) {
            while(p[i].x[d] < 0) p[i].x[d] += BoxSize;
            while(p[i].x[d] >= BoxSize) p[i].x[d] -= BoxSize;
        }
    }
}

static int domain_decompose0(void * P, size_t elsize, int np_local, int np_allocated, 
        void * scratch, size_t scratch_bytes);

void domain_decompose(Particles* particles, void * scratch, size_t scratch_bytes) {
    void * P = particles->p;
    size_t elsize = sizeof(Particle);
    particles->np_local = 
        domain_decompose0(P, elsize, particles->np_local, particles->np_allocated, scratch, scratch_bytes);
}

void domain_decompose_min(Snapshot * particles, void * scratch, size_t scratch_bytes) {
    void * P = particles->p;
    size_t elsize = sizeof(ParticleMinimum);
    particles->np_local = 
        domain_decompose0(P, elsize, particles->np_local, particles->np_allocated, scratch, scratch_bytes);
}

int domain_decompose0(void * P, size_t elsize, int np_local, int np_allocated, 
        void * scratch, size_t scratch_bytes)
{
    char * Pc = P;

    {
    size_t total1 = np_local;
    size_t total = 0;
    MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Total number of particles = %td\n", total);
    total1 = 0;
    for(int i = 0; i < np_local; i ++) {
        total1 += (((ParticleMinimum*) &Pc[i * elsize])->id);
    }
    MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Sum of particle ID = %td\n", total);
    }

    qsort_r_with_tmp(P, np_local, elsize, cmp_par_task, NULL, scratch, scratch_bytes);

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
    for(int i = 0; i < np_local; i ++) {
        int task = par_node(((ParticleMinimum*) &Pc[i * elsize])->x);
        SendCount[task] ++;
    }
    char * sendbuf = (char*) &Pc[SendCount[ThisTask] *elsize];
    char * recvbuf = scratch;

    MPI_Datatype MPI_PAR;
    MPI_Type_contiguous(elsize, MPI_BYTE, &MPI_PAR);
    MPI_Type_commit(&MPI_PAR);

    /* chop the p buffer to local and send buf */
    np_local = SendCount[ThisTask];
    /* local particles are not sent */
    SendCount[ThisTask] = 0;

    MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
    }
    /* prepare the recv buf from scratch */
    size_t recv_total = RecvDispl[NTask - 1] + RecvCount[NTask - 1];
    if(recv_total * elsize > scratch_bytes) {
        msg_abort(2314, "Not enough memory in scratch for this exchange."
                "Need %td bytes. Provided %td bytes\n",
                recv_total, scratch_bytes);
    }
    
    MPI_Alltoallv(sendbuf, SendCount, SendDispl, MPI_PAR, 
            recvbuf, RecvCount, RecvDispl, MPI_PAR, MPI_COMM_WORLD);
    
    if(recv_total + np_local > np_allocated) {
        msg_abort(2314, "Not enough memory in particles for this exchange."
                "Need %d particles. Allocated %d particles\n",
                recv_total + np_local, np_allocated);

    }
    /* concatenate received particles with the local particles */
    memcpy(sendbuf, recvbuf, elsize * recv_total);
    np_local += recv_total;

    MPI_Type_free(&MPI_PAR);
    for(int i = 0; i < np_local; i ++) {
        if(par_node(((ParticleMinimum*) &Pc[i * elsize])->x) != ThisTask) {
            msg_abort(2314, "Some particles are not in the correct domain after exchange.");
        }
    }
    {
    size_t total1 = np_local;
    size_t total = 0;
    MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Total number of particles = %td\n", total);
    total1 = 0;
    for(int i = 0; i < np_local; i ++) {
        total1 += (((ParticleMinimum*) &Pc[i * elsize])->id);
    }
    MPI_Allreduce(&total1, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    msg_printf(verbose, "Sum of particle ID = %td\n", total);
    }
    return np_local;
}
