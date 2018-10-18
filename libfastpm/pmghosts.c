#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"
#include "pmghosts.h"

#ifdef ENABLE_VALGRIND
#include </usr/include/valgrind/memcheck.h>
#endif

typedef void (*pm_iter_ghosts_func)(PM * pm, PMGhostData * ppd, void * userdata);

void pm_ghosts_free(PMGhostData * pgd) {
    if(pgd->p) {
        fastpm_store_destroy(pgd->p);
        free(pgd->p);
    }

    fastpm_memory_free(pgd->pm->mem, pgd->ighost_to_ipar);

    free(pgd->Nsend);
    free(pgd->Osend);
    free(pgd->Nrecv);
    free(pgd->Orecv);
    free(pgd);
}

static void
pm_iter_ghosts(PM * pm, PMGhostData * pgd, 
    pm_iter_ghosts_func iter_func, void * userdata)
{

    ptrdiff_t i;
    for (i = 0; i < pgd->source->np; i ++) {
        PMGhostData localppd = *pgd;
        double pos[3];
        int rank;
        fastpm_store_get_position(pgd->source, i, pos);
        int d;

        /* how far the window expands. */
        int left[3];
        int right[3];
        for(d = 0; d < 3; d ++) {
            /* this condition is not tightest for CIC painting, because
             * a particle touches a cell doesn't mean cic touches the left edge
             * of the cell.
             * */
            left[d] = floor(pos[d] * pm->InvCellSize[d] + pgd->Below[d]);
            right[d] = floor(pos[d] * pm->InvCellSize[d] + pgd->Above[d]);
        }

        /* probe neighbours */
        int j[3];
        int ranks[1000];
        int used = 0;
        localppd.ipar = i;
        /* no need to run the z loop because the decomposition is in xy */
        for(j[2] = left[2]; j[2] <= right[2]; j[2] ++)
        for(j[0] = left[0]; j[0] <= right[0]; j[0] ++)
        for(j[1] = left[1]; j[1] <= right[1]; j[1] ++)
        {
            rank = pm_ipos_to_rank(pm, j);
            if(LIKELY(rank == pm->ThisTask))  continue;
            int ptr;
            for(ptr = 0; ptr < used; ptr++) {
                if(rank == ranks[ptr]) break;
            } 
            if(UNLIKELY(ptr == used)) {
                ranks[used++] = rank;
                localppd.rank = rank;
                localppd.reason = j;
                iter_func(pm, &localppd, userdata);
            }
        }
    }
}

static void
count_ghosts(PM * pm, PMGhostData * pgd, void * userdata)
{
#pragma omp atomic
    pgd->Nsend[pgd->rank] ++;
}

struct build_ghost_buffer_data {
    enum FastPMPackFields attributes;
    size_t elsize;
};

static void
build_ghost_buffer(PM * pm, PMGhostData * pgd, void * userdata)
{
    struct build_ghost_buffer_data * data = userdata;


    int ighost;
    int offset; 

#pragma omp atomic capture
    offset = pgd->Nsend[pgd->rank] ++;

    ighost = pgd->Osend[pgd->rank] + offset;

    fastpm_store_pack(pgd->source, pgd->ipar,
        (char*) pgd->send_buffer + ighost * data->elsize, data->attributes);

    pgd->ighost_to_ipar[ighost] = pgd->ipar;
}

/* create ghosts that can hold 'attributes';
 * use pm_ghosts_send to send subsets of `attributes`;
 * */
PMGhostData *
pm_ghosts_create(PM * pm, FastPMStore * p,
    enum FastPMPackFields attributes)
{
    return pm_ghosts_create_full(pm, p, attributes, pm->Below, pm->Above);

}

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        enum FastPMPackFields attributes,
        double below[],
        double above[]
        )
{
    PMGhostData * pgd = malloc(sizeof(pgd[0]));

    pgd->pm = pm;
    pgd->source = p;

    int d;
    for(d = 0; d < 3; d++) {
        pgd->Below[d] = below[d] * pm->InvCellSize[d];
        pgd->Above[d] = above[d] * pm->InvCellSize[d];
    }

    pgd->ighost_to_ipar = NULL;

    pgd->Nsend = calloc(pm->NTask, sizeof(int));
    pgd->Osend = calloc(pm->NTask, sizeof(int));
    pgd->Nrecv = calloc(pm->NTask, sizeof(int));
    pgd->Orecv = calloc(pm->NTask, sizeof(int));

    size_t Nsend;
    size_t Nrecv;

    memset(pgd->Nsend, 0, sizeof(pgd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, pgd, count_ghosts, NULL);

    Nsend = cumsum(pgd->Osend, pgd->Nsend, pm->NTask);

    MPI_Alltoall(pgd->Nsend, 1, MPI_INT, pgd->Nrecv, 1, MPI_INT, pm->Comm2D);

    Nrecv = cumsum(pgd->Orecv, pgd->Nrecv, pm->NTask);

    pgd->ighost_to_ipar = fastpm_memory_alloc(pm->mem, "Ghost2Par", Nsend * sizeof(int), FASTPM_MEMORY_HEAP);

    pgd->p = malloc(sizeof(pgd->p[0]));
    fastpm_store_init(pgd->p, Nrecv, attributes, FASTPM_MEMORY_HEAP);

    return pgd;
}

void
pm_ghosts_send(PMGhostData * pgd, enum FastPMPackFields attributes)
{
    PM * pm = pgd->pm;
    ptrdiff_t i;
    size_t Nsend;
    size_t Nrecv;

    Nsend = cumsum(pgd->Osend, pgd->Nsend, pm->NTask);
    Nrecv = cumsum(pgd->Orecv, pgd->Nrecv, pm->NTask);

    size_t elsize = fastpm_store_pack(pgd->source, 0, NULL, attributes);

    pgd->send_buffer = fastpm_memory_alloc(pm->mem, "SendBuf", Nsend * elsize, FASTPM_MEMORY_STACK);
    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, "RecvBuf", Nrecv * elsize, FASTPM_MEMORY_STACK);

    /* build buffer */
    memset(pgd->Nsend, 0, sizeof(pgd->Nsend[0]) * pm->NTask);

    struct build_ghost_buffer_data data = {attributes, elsize};
    pm_iter_ghosts(pm, pgd, build_ghost_buffer, &data);

    /* exchange */

    pgd->p->np = Nrecv;

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(pgd->send_buffer, pgd->Nsend, pgd->Osend, GHOST_TYPE,
                  pgd->recv_buffer, pgd->Nrecv, pgd->Orecv, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

#pragma omp parallel for
    for(i = 0; i < Nrecv; i ++) {
        fastpm_store_unpack(pgd->p, i,
                (char*) pgd->recv_buffer + i * elsize,
                        attributes);
    }
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
    fastpm_memory_free(pm->mem, pgd->send_buffer);
}

void
pm_ghosts_reduce(PMGhostData * pgd, FastPMFieldDescr field)
{

    int ci;
    for(ci = 0; ci < 32; ci ++) {
        if(field.attribute == pgd->p->column_info[ci].attribute) break;
    }

    PM * pm = pgd->pm;

    size_t Nsend = cumsum(NULL, pgd->Nsend, pm->NTask);
    size_t Nrecv = cumsum(NULL, pgd->Nrecv, pm->NTask);
    ptrdiff_t i;

    size_t elsize;

    elsize = pgd->p->column_info[ci].elsize;

    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, "RecvBuf", Nrecv * elsize, FASTPM_MEMORY_STACK);
    pgd->send_buffer = fastpm_memory_alloc(pm->mem, "SendBuf", Nsend * elsize, FASTPM_MEMORY_STACK);

#pragma omp parallel for
    for(i = 0; i < pgd->p->np; i ++) {
        pgd->p->column_info[ci].pack_member(pgd->p, i, ci, field.memb,
            (char*) pgd->recv_buffer + i * elsize);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(pgd->recv_buffer, pgd->Nrecv, pgd->Orecv, GHOST_TYPE,
                  pgd->send_buffer, pgd->Nsend, pgd->Osend, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

    /* now reduce the attributes. */
    int ighost;

    /* this loop is not parallel because multiple ghosts can be for the same ipar,
     * in which case we have a race condition.
     * we can fix this by carefully working with ipar (it should / could be made sorted)
     * but unlikly worth the effort.
     * */
    for(ighost = 0; ighost < Nsend; ighost ++) {
        pgd->p->column_info[ci].reduce_member(pgd->source,
            pgd->ighost_to_ipar[ighost], ci, field.memb,
            pgd->send_buffer + ighost * elsize);
    }
    fastpm_memory_free(pm->mem, pgd->send_buffer);
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
}
