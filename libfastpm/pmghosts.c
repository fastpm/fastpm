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

static void
build_ghost_buffer(PM * pm, PMGhostData * pgd, void * userdata)
{
    FastPMPackingPlan * plan = userdata;


    int ighost;
    int offset; 

#pragma omp atomic capture
    offset = pgd->Nsend[pgd->rank] ++;

    ighost = pgd->Osend[pgd->rank] + offset;

    fastpm_packing_plan_pack(plan, pgd->source, pgd->ipar, 
                (char*) pgd->send_buffer + ighost * plan->elsize);

    pgd->ighost_to_ipar[ighost] = pgd->ipar;
}

/* create ghosts that can hold 'attributes';
 * use pm_ghosts_send to send subsets of `attributes`;
 * */
PMGhostData *
pm_ghosts_create(PM * pm, FastPMStore * p,
    FastPMColumnTags attributes, int support)
{
    /* The support of CIC is 2. We do not use
     * -1.0000 * pm->CellSize[d] here
     * because even though the kernel touches -1 * cellsize,
     * we do not paint on the lower edge.
     * */
    double Below[3]; /* in grid integer units */
    double Above[3];

    int d;
    for(d = 0; d < 3; d ++) {
        Below[d] = - (support * 0.5 - 1);
        Above[d] = (support * 0.5    );
    }
    return pm_ghosts_create_full(pm, p, attributes, Below, Above);

}

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        FastPMColumnTags attributes,
        double below[],
        double above[]
        )
{
    PMGhostData * pgd = malloc(sizeof(pgd[0]));

    pgd->pm = pm;
    pgd->source = p;

    int d;
    for(d = 0; d < 3; d++) {
        pgd->Below[d] = below[d];
        pgd->Above[d] = above[d];
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

    double nmin, nmax, nmean, nstd;
    MPIU_stats(pm->Comm2D, Nsend, "<->s", &nmin, &nmean, &nmax, &nstd);
    fastpm_info("Sending ghosts: min = %g max = %g mean = %g std = %g\n",
        nmin, nmax, nmean, nstd);

    MPIU_stats(pm->Comm2D, Nrecv, "<->s", &nmin, &nmean, &nmax, &nstd);
    fastpm_info("Receiving ghosts: min = %g max = %g mean = %g std = %g\n",
        nmin, nmax, nmean, nstd);

    pgd->ighost_to_ipar = fastpm_memory_alloc(pm->mem, "Ghost2Par", Nsend * sizeof(int), FASTPM_MEMORY_HEAP);

    pgd->p = malloc(sizeof(pgd->p[0]));
    fastpm_store_init(pgd->p, pgd->source->name, Nrecv, attributes, FASTPM_MEMORY_HEAP);
    memcpy(&pgd->p->meta, &pgd->source->meta, sizeof(pgd->source->meta));

    return pgd;
}

void
pm_ghosts_has_ghosts(PMGhostData * pgd, uint8_t * has_ghosts)
{
    size_t Nsend = cumsum(NULL, pgd->Nsend, pgd->pm->NTask);

    ptrdiff_t i;
    for(i = 0; i < pgd->source->np; i ++) {
        has_ghosts[i] = 0;
    }
    for(i = 0; i < Nsend; i ++) {
        has_ghosts[pgd->ighost_to_ipar[i]] = 1;
    }
}

void
pm_ghosts_send(PMGhostData * pgd, FastPMColumnTags attributes)
{
    PM * pm = pgd->pm;
    ptrdiff_t i;
    size_t Nsend;
    size_t Nrecv;

    Nsend = cumsum(pgd->Osend, pgd->Nsend, pm->NTask);
    Nrecv = cumsum(pgd->Orecv, pgd->Nrecv, pm->NTask);

    FastPMPackingPlan plan[1];
    fastpm_packing_plan_init(plan, pgd->p, attributes);

    pgd->send_buffer = fastpm_memory_alloc(pm->mem, "SendBuf", Nsend * plan->elsize, FASTPM_MEMORY_STACK);
    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, "RecvBuf", Nrecv * plan->elsize, FASTPM_MEMORY_STACK);

    /* build buffer */
    memset(pgd->Nsend, 0, sizeof(pgd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, pgd, build_ghost_buffer, plan);

    /* exchange */

    pgd->p->np = Nrecv;

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(plan->elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(pgd->send_buffer, pgd->Nsend, pgd->Osend, GHOST_TYPE,
                  pgd->recv_buffer, pgd->Nrecv, pgd->Orecv, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

#pragma omp parallel for
    for(i = 0; i < Nrecv; i ++) {
        fastpm_packing_plan_unpack(plan,
                pgd->p, i,
                (char*) pgd->recv_buffer + i * plan->elsize);
    }
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
    fastpm_memory_free(pm->mem, pgd->send_buffer);
}

void
pm_ghosts_reduce(PMGhostData * pgd, FastPMColumnTags attribute,
    reduce_func reduce,
    void * userdata
)
{

    int ci = fastpm_store_find_column_id(pgd->p, attribute);

    PM * pm = pgd->pm;

    size_t Nsend = cumsum(NULL, pgd->Nsend, pm->NTask);
    size_t Nrecv = cumsum(NULL, pgd->Nrecv, pm->NTask);
    ptrdiff_t i;

    size_t elsize;

    elsize = pgd->p->_column_info[ci].elsize;

    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, "RecvBuf", Nrecv * elsize, FASTPM_MEMORY_STACK);
    pgd->send_buffer = fastpm_memory_alloc(pm->mem, "SendBuf", Nsend * elsize, FASTPM_MEMORY_STACK);

#pragma omp parallel for
    for(i = 0; i < pgd->p->np; i ++) {
        pgd->p->_column_info[ci].pack(pgd->p, i, ci,
            (char*) pgd->recv_buffer + i * elsize);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(pgd->recv_buffer, pgd->Nrecv, pgd->Orecv, GHOST_TYPE,
                  pgd->send_buffer, pgd->Nsend, pgd->Osend, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

    FastPMStore q[1];
    fastpm_store_init(q, pgd->p->name, Nsend, attribute, FASTPM_MEMORY_HEAP);

    /* now reduce the attributes. */
    int ighost;

    /* this loop is not parallel because multiple ghosts can be for the same ipar,
     * in which case we have a race condition.
     * we can fix this by carefully working with ipar (it should / could be made sorted)
     * but unlikly worth the effort.
     * */
    for(ighost = 0; ighost < Nsend; ighost ++) {
        pgd->p->_column_info[ci].unpack(q, ighost, ci,
            (char*) pgd->send_buffer + ighost * elsize);
    }

    for(ighost = 0; ighost < Nsend; ighost ++) {
        reduce(q, ighost, pgd->source, pgd->ighost_to_ipar[ighost],
                ci, userdata);
    }

    fastpm_store_destroy(q);
    fastpm_memory_free(pm->mem, pgd->send_buffer);
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
}
