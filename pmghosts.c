#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include "pmpfft.h"
#include "msg.h"

void pm_destroy_ghosts(PMGhostData * ppd) {
    ppd->pm->iface.free(ppd->ighost_to_ipar);
    free(ppd->Nsend);
    free(ppd->Osend);
    free(ppd->Nrecv);
    free(ppd->Orecv);
}

static void pm_iter_ghosts(PM * pm, PMGhostData * ppd, 
    pm_iter_ghosts_func iter_func) {

    ptrdiff_t i;
    for (i = 0; i < ppd->np; i ++) {
        PMGhostData localppd = *ppd;
        double pos[3];
        int rank;
        ppd->get_position(ppd->pdata, i, pos);
        int d;
        int ipos[3];
        for(d = 0; d < 3; d ++) {
            ipos[d] = floor(pos[d] * pm->InvCellSize[d]);
        }

        /* probe neighbours */
        ptrdiff_t j[3];
        int ranks[1000];
        int used = 0;
        localppd.ipar = i;
        for(j[2] = pm->Below[2]; j[2] <= pm->Above[2]; j[2] ++)
        for(j[0] = pm->Below[0]; j[0] <= pm->Above[0]; j[0] ++)
        for(j[1] = pm->Below[1]; j[1] <= pm->Above[1]; j[1] ++) {
            int npos[3];
            int d;
            for(d = 0; d < 3; d ++) {
                npos[d] = ipos[d] + j[d];
            }
            rank = pm_ipos_to_rank(pm, npos);
            if(LIKELY(rank == pm->ThisTask))  continue;
            int ptr;
            for(ptr = 0; ptr < used; ptr++) {
                if(rank == ranks[ptr]) break;
            } 
            if(UNLIKELY(ptr == used)) {
                ranks[used++] = rank;
                localppd.rank = rank;
                localppd.reason = j;
                iter_func(pm, &localppd);
            } 
        }
    }
}

static void count_ghosts(PM * pm, PMGhostData * ppd) {
//#pragma omp atomic
    ppd->Nsend[ppd->rank] ++;
}

static void build_ghost_buffer(PM * pm, PMGhostData * ppd) {
    double pos[3];
    ppd->get_position(ppd->pdata, ppd->ipar, pos);

    int ighost;
    int offset; 

//#pragma omp atomic capture
    offset = ppd->Nsend[ppd->rank] ++;

    ighost = ppd->Osend[ppd->rank] + offset;

    pm->iface.pack(ppd->pdata, ppd->ipar, 
        (char*) ppd->send_buffer + ighost * ppd->elsize, ppd->attributes);

    ppd->ighost_to_ipar[ighost] = ppd->ipar;
#if 0
    msg_aprintf(debug, "Making a ghost for particle %td (%g %g %g) to rank %d for %td %td %td\n", 
            ppd->ipar, 
            pos[0], pos[1], pos[2],
            ppd->rank, ppd->reason[0], ppd->reason[1], ppd->reason[2]);

    msg_aprintf(debug, "Connecting Ghost %d to Particle %td\n", ighost, ppd->ipar);
#endif
}

void pm_append_ghosts(PMGhostData * ppd) {
    PM * pm = ppd->pm;
    ptrdiff_t i;
    size_t Nsend;
    size_t Nrecv;
    size_t elsize = pm->iface.pack(ppd->pdata, 0, NULL, ppd->attributes);

    ppd->Nsend = calloc(pm->NTask, sizeof(int));
    ppd->Osend = calloc(pm->NTask, sizeof(int));
    ppd->Nrecv = calloc(pm->NTask, sizeof(int));
    ppd->Orecv = calloc(pm->NTask, sizeof(int));

    ppd->elsize = elsize;

    memset(ppd->Nsend, 0, sizeof(ppd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, ppd, count_ghosts);

    Nsend = cumsum(ppd->Osend, ppd->Nsend, pm->NTask);

    MPI_Alltoall(ppd->Nsend, 1, MPI_INT, ppd->Nrecv, 1, MPI_INT, pm->Comm2D);

    Nrecv = cumsum(ppd->Orecv, ppd->Nrecv, pm->NTask);
    

    ppd->ighost_to_ipar = pm->iface.malloc(Nsend * sizeof(int));
    ppd->send_buffer = pm->iface.malloc(Nsend * ppd->elsize);
    ppd->recv_buffer = pm->iface.malloc(Nrecv * ppd->elsize);

    memset(ppd->Nsend, 0, sizeof(ppd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, ppd, build_ghost_buffer);

    /* exchange */

    ppd->nghosts = Nrecv;

    if(Nrecv + ppd->np > ppd->np_upper) {
        msg_abort(-1, "Too many ghosts; asking for %td, space for %td\n", Nrecv, ppd->np_upper - ppd->np);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(ppd->elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(ppd->send_buffer, ppd->Nsend, ppd->Osend, GHOST_TYPE,
                  ppd->recv_buffer, ppd->Nrecv, ppd->Orecv, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

#pragma omp parallel for
    for(i = 0; i < Nrecv; i ++) {
        pm->iface.unpack(ppd->pdata, ppd->np + i, 
                (char*) ppd->recv_buffer + i * ppd->elsize, 
                        ppd->attributes);
    }
    pm->iface.free(ppd->recv_buffer);
    pm->iface.free(ppd->send_buffer);
}

void pm_reduce_ghosts(PMGhostData * ppd, int attributes) {
    PM * pm = ppd->pm;
    size_t Nsend = cumsum(NULL, ppd->Nsend, pm->NTask);
    size_t Nrecv = cumsum(NULL, ppd->Nrecv, pm->NTask);
    ptrdiff_t i;

    ppd->elsize = pm->iface.pack(ppd->pdata, 0, NULL, attributes);
    ppd->recv_buffer = pm->iface.malloc(Nrecv * ppd->elsize);
    ppd->send_buffer = pm->iface.malloc(Nsend * ppd->elsize);
    ppd->ReductionAttributes = attributes;

#pragma omp parallel for
    for(i = 0; i < ppd->nghosts; i ++) {
        pm->iface.pack(ppd->pdata, i + ppd->np, 
            (char*) ppd->recv_buffer + i * ppd->elsize, 
            ppd->ReductionAttributes);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(ppd->elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(ppd->recv_buffer, ppd->Nrecv, ppd->Orecv, GHOST_TYPE,
                  ppd->send_buffer, ppd->Nsend, ppd->Osend, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

    /* now reduce the attributes. */
    int ighost;
#pragma omp parallel for
    for(ighost = 0; ighost < Nsend; ighost ++) {
        pm->iface.reduce(ppd->pdata, ppd->ighost_to_ipar[ighost], 
            (char*) ppd->send_buffer + ighost * ppd->elsize, 
            ppd->ReductionAttributes);
    }
    pm->iface.free(ppd->send_buffer);
    pm->iface.free(ppd->recv_buffer);
}
