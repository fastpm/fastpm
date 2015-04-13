#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <assert.h>
#include <math.h>
#include "msg.h"
#include "move.h"
#include "comm.h"
#include "particle.h"
#include <fftw3-mpi.h>

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

void domain_init(int Nmesh_, double BoxSize) {
    Nmesh = Nmesh_;
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
static int par_node(const ParticleMinimum * p) {
    int task = 0;
    for(int d = 0; d < 2; d ++) {
        int ind = floor(p->x[d] * MeshPerBox);
        while(ind < 0) ind += Nmesh;
        while(ind >= Nmesh) ind -= Nmesh;
        task += TaskStrides[d] * MeshToTask[d][ind];
    }
    return task;
}

/* this will sort my particles to the top, 
 * then the other particles near the end by order */
static int cmp_par_task(const void * c1, const void * c2) {
    int task1 = par_node((ParticleMinimum*) c1);
    int task2 = par_node((ParticleMinimum*) c2);
    if (task1 == ThisTask) task1 = -1;
    if (task2 == ThisTask) task2 = -1;
    return (task1 > task2) - (task2 > task1);
}

void domain_decompose(Particles* const particles, void * scratch, size_t scratch_bytes)
{
    size_t elsize = sizeof(Particle);
    qsort(particles->p, particles->np_local, elsize, cmp_par_task);

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
    for(int i = 0; i < particles->np_local; i ++) {
        int task = par_node((ParticleMinimum*) &particles->p[i]);
        SendCount[task] ++;
    }
    char * sendbuf = (char*) &particles->p[SendCount[ThisTask]];
    char * recvbuf = scratch;

    MPI_Datatype MPI_PAR;
    MPI_Type_contiguous(elsize, MPI_BYTE, &MPI_PAR);
    MPI_Type_commit(&MPI_PAR);

    /* chop the p buffer to local and send buf */
    particles->np_local = SendCount[ThisTask];
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
    
    if(recv_total + particles->np_local > particles->np_allocated) {
        msg_abort(2314, "Not enough memory in particles for this exchange."
                "Need %d particles. Allocated %d particles\n",
                recv_total + particles->np_local, particles->np_allocated);

    }
    /* concatenate received particles with the local particles */
    memcpy(sendbuf, recvbuf, elsize * recv_total);
    particles->np_local += recv_total;

    MPI_Type_free(&MPI_PAR);
    for(int i = 0; i < particles->np_local; i ++) {
        if(par_node((ParticleMinimum*) &particles->p[i]) != ThisTask) {
            msg_abort(2314, "Some particles are not in the correct domain after exchange.");
        }
    }
}
