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
    for (i = 0; i < pgd->np; i ++) {
        PMGhostData localppd = *pgd;
        double pos[3];
        int rank;
        pgd->get_position(pgd->p, i, pos);
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
    enum FastPMPackFields * attributes = userdata;

    FastPMStore * p = pgd->p;
    double pos[3];
    pgd->get_position(pgd->p, pgd->ipar, pos);

    int ighost;
    int offset; 

#pragma omp atomic capture
    offset = pgd->Nsend[pgd->rank] ++;

    ighost = pgd->Osend[pgd->rank] + offset;

    p->pack(p, pgd->ipar,
        (char*) pgd->send_buffer + ighost * pgd->elsize, *attributes);

    pgd->ighost_to_ipar[ighost] = pgd->ipar;
}

PMGhostData *
pm_ghosts_create(PM * pm, FastPMStore *p,
    enum FastPMPackFields attributes,
    fastpm_posfunc get_position)
{
    return
    pm_ghosts_create_full(pm, p, attributes,
            get_position, pm->Below, pm->Above);
}

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        enum FastPMPackFields attributes,
        fastpm_posfunc get_position,
        double below[],
        double above[]
        )
{
    PMGhostData * pgd = malloc(sizeof(pgd[0]));

    pgd->pm = pm;
    pgd->p = p;
    pgd->np = p->np;
    pgd->np_upper = p->np_upper;

    if(get_position == NULL)
        pgd->get_position = p->get_position;
    else
        pgd->get_position = get_position;
    pgd->nghosts = 0;

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

    pm_ghosts_send(pgd, attributes);

    return pgd;
}

void
pm_ghosts_send(PMGhostData * pgd, enum FastPMPackFields attributes)
{
    PM * pm = pgd->pm;
    FastPMStore * p = pgd->p;

    ptrdiff_t i;
    size_t Nsend;
    size_t Nrecv;
    size_t elsize = p->pack(pgd->p, 0, NULL, attributes);

    pgd->elsize = elsize;

    memset(pgd->Nsend, 0, sizeof(pgd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, pgd, count_ghosts, NULL);

    Nsend = cumsum(pgd->Osend, pgd->Nsend, pm->NTask);

    MPI_Alltoall(pgd->Nsend, 1, MPI_INT, pgd->Nrecv, 1, MPI_INT, pm->Comm2D);

    Nrecv = cumsum(pgd->Orecv, pgd->Nrecv, pm->NTask);

    if(pgd->ighost_to_ipar == NULL)
        pgd->ighost_to_ipar = fastpm_memory_alloc(pm->mem, Nsend * sizeof(int), FASTPM_MEMORY_HEAP);

    pgd->send_buffer = fastpm_memory_alloc(pm->mem, Nsend * pgd->elsize, FASTPM_MEMORY_HEAP);
    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, Nrecv * pgd->elsize, FASTPM_MEMORY_HEAP);

    memset(pgd->Nsend, 0, sizeof(pgd->Nsend[0]) * pm->NTask);

    pm_iter_ghosts(pm, pgd, build_ghost_buffer, &attributes);

    /* exchange */

    pgd->nghosts = Nrecv;

    if(Nrecv + pgd->np > pgd->np_upper) {
        fastpm_raise(-1, "Too many ghosts; asking for %td, space for %td\n", Nrecv, pgd->np_upper - pgd->np);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(pgd->elsize, MPI_BYTE, &GHOST_TYPE);
    MPI_Type_commit(&GHOST_TYPE);
    MPI_Alltoallv_sparse(pgd->send_buffer, pgd->Nsend, pgd->Osend, GHOST_TYPE,
                  pgd->recv_buffer, pgd->Nrecv, pgd->Orecv, GHOST_TYPE,
                    pm->Comm2D);
    MPI_Type_free(&GHOST_TYPE);

#pragma omp parallel for
    for(i = 0; i < Nrecv; i ++) {
        p->unpack(pgd->p, pgd->np + i,
                (char*) pgd->recv_buffer + i * pgd->elsize,
                        attributes);
    }
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
    fastpm_memory_free(pm->mem, pgd->send_buffer);
}

static
void
_reduce_any_field(PMGhostData * pgd,
    enum FastPMPackFields attributes,
    ptrdiff_t index, void * buf,
    void * userdata)
{
    pgd->p->reduce(pgd->p, index, buf, attributes);
}

void
pm_ghosts_reduce(PMGhostData * pgd, enum FastPMPackFields attributes)
{
    pm_ghosts_reduce_any(pgd, attributes,
        (pm_ghosts_reduce_func) _reduce_any_field,
        NULL);
}

void
pm_ghosts_reduce_any(PMGhostData * pgd, 
        enum FastPMPackFields attributes,
        pm_ghosts_reduce_func func,
        void * userdata)
{
    PM * pm = pgd->pm;
    FastPMStore * p = pgd->p;

    size_t Nsend = cumsum(NULL, pgd->Nsend, pm->NTask);
    size_t Nrecv = cumsum(NULL, pgd->Nrecv, pm->NTask);
    ptrdiff_t i;

    pgd->elsize = p->pack(pgd->p, 0, NULL, attributes);
    pgd->recv_buffer = fastpm_memory_alloc(pm->mem, Nrecv * pgd->elsize, FASTPM_MEMORY_HEAP);
    pgd->send_buffer = fastpm_memory_alloc(pm->mem, Nsend * pgd->elsize, FASTPM_MEMORY_HEAP);

#pragma omp parallel for
    for(i = 0; i < pgd->nghosts; i ++) {
        p->pack(p, i + pgd->np,
            (char*) pgd->recv_buffer + i * pgd->elsize, 
            attributes);
    }

    MPI_Datatype GHOST_TYPE;
    MPI_Type_contiguous(pgd->elsize, MPI_BYTE, &GHOST_TYPE);
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
        func(pgd, attributes, 
            pgd->ighost_to_ipar[ighost],
            pgd->send_buffer + ighost * pgd->elsize,
            userdata);
    }
    fastpm_memory_free(pm->mem, pgd->send_buffer);
    fastpm_memory_free(pm->mem, pgd->recv_buffer);
}
