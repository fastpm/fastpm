#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

typedef struct FastPMDomain {
    FastPMMemory * mem;
    PM * pm;
    ptrdiff_t * masterindex;
    MPI_Comm comm;
    size_t oldsize;
    size_t newsize;

} FastPMDomain;

void 
fastpm_domain_init(FastPMDomain * self, PM * pm, FastPMColumn * pos, MPI_Comm comm)
{
    self->mem = pos->mem;
    self->pm = pm;
    self->comm = comm;
    self->masterindex = fastpm_memory_alloc(self->mem, pos->size * sizeof(self->masterindex[0]), FASTPM_MEMORY_HEAP);

    /* decompose, and create ghost data */
    double x[3];
    ptrdiff_t i;

    int ThisTask;
    int NTask;
    MPI_Comm_rank(self->comm, &ThisTask);
    MPI_Comm_size(self->comm, &NTask);

    long long * offset = malloc(sizeof(long long) * (NTask + 1));
    long long * counts = malloc(sizeof(long long) * (NTask));

    offset[ThisTask + 1] = pos->size;
    offset[0] = 0;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_LONG_LONG, offset + 1, 1, MPI_LONG_LONG, self->comm);

    for(i = 1; i <= NTask; i ++) {
        offset[i] += offset[i - 1];
    }
    for(i = 0; i < NTask; i ++) {
        counts[i] = 0;
    }
    for(i = 0; i < pos->size; i ++) {
        fastpm_column_get_double(pos, i, x);
        int ipos[3];
        int masterrank = pm_pos_to_rank(self->pm, x, ipos);
        self->masterindex[i] = masterrank * offset[NTask] + i + offset[ThisTask];

        counts[masterrank] ++;
    }
    MPI_Allreduce(MPI_IN_PLACE, counts, NTask, MPI_LONG_LONG, MPI_SUM, self->comm);

    self->oldsize = pos->size;
    self->newsize = counts[ThisTask];

    free(offset);
    free(counts);
}

void
fastpm_domain_decompose(FastPMDomain * self, FastPMColumn * columns[])
{
    FastPMColumnSet columnset[1];
    fastpm_columnset_init(columnset, columns);

    fastpm_column_parallel_permute((FastPMColumn*)columnset, self->masterindex, self->newsize, self->comm);

    fastpm_column_destroy((FastPMColumn*) columnset);
}

void
fastpm_domain_destroy(FastPMDomain * self)
{
    fastpm_memory_free(self->mem, self->masterindex);
}

struct FastPMGhostPair {
    int ghostrank;
    int i;
};
static int
_ghost_data_sort_by_rank(const void * p1, const void * p2)
{
    const int * i1 = p1;
    const int * i2 = p2;
    return *i1 - *i2;
}

typedef struct FastPMGhosts {
    FastPMMemory * mem;
    PM * pm;
    double margin;
    struct FastPMGhostPair * ghostindex;
    MPI_Comm comm;

    size_t nsend;
    size_t nrecv;
    int * sendcounts;
    int * senddispls;
    int * recvcounts;
    int * recvdispls;
} FastPMGhosts;

static int
 _fastpm_ghosts_probe(FastPMGhosts * self, double x[3], int masterrank, int * ranks);

void
fastpm_ghosts_init(FastPMGhosts * self, double margin, PM * pm, FastPMColumn * pos, MPI_Comm comm)
{
    self->mem = pos->mem;
    self->comm = comm;
    self->pm = pm;
    self->margin = margin;
    int ThisTask;
    int NTask;
    MPI_Comm_rank(self->comm, &ThisTask);
    MPI_Comm_size(self->comm, &NTask);

    self->sendcounts = malloc(sizeof(int) * NTask);
    self->recvcounts = malloc(sizeof(int) * NTask);
    self->senddispls = malloc(sizeof(int) * NTask);
    self->recvdispls = malloc(sizeof(int) * NTask);

    int * ranks = (int*) malloc(sizeof(int) * NTask);

    size_t maxghosts = pos->size;
    while(1) {
        self->ghostindex = fastpm_memory_alloc(self->mem, sizeof(self->ghostindex[0]) * maxghosts, FASTPM_MEMORY_HEAP);
        ptrdiff_t i;
        double x[3];
        ptrdiff_t gi = 0;

        for(i = 0; i < pos->size; i ++) {
            fastpm_column_get_double(pos, i, x);
            int nranks = _fastpm_ghosts_probe(self, x, ThisTask, ranks);
            int j;
            for(j = 0; j < nranks; j ++) {
                self->ghostindex[gi].ghostrank = ranks[j];
                self->ghostindex[gi].i = i;
                gi ++;
                if(gi >= maxghosts) {
                    goto bail;
                }
            }
        }
        self->nsend = gi;
        break;
    bail:
        maxghosts = maxghosts + ((maxghosts + 1) >> 1);
        fastpm_memory_free(self->mem, self->ghostindex);
    }

    qsort(self->ghostindex, self->nsend, sizeof(self->ghostindex[0]), _ghost_data_sort_by_rank);

    long long * sendcounts = malloc(sizeof(long long) * (NTask));
    long long * recvcounts = malloc(sizeof(long long) * (NTask));
    long long * senddispls = malloc(sizeof(long long) * (NTask));
    long long * recvdispls = malloc(sizeof(long long) * (NTask));

    ptrdiff_t i;

    for(i = 0; i < self->nsend; i ++) {
        sendcounts[self->ghostindex[i].ghostrank] ++;
    }
    MPI_Alltoall(sendcounts, 1, MPI_LONG_LONG, recvcounts, 1, MPI_LONG_LONG, self->comm);

    recvdispls[0] = 0;
    senddispls[0] = 0;
    for(i = 1; i < NTask; i ++) {
        senddispls[i] = senddispls[i - 1] + sendcounts[i];
        recvdispls[i] = recvdispls[i - 1] + recvcounts[i];
    }

    for(i = 0; i < NTask; i ++) {
        self->sendcounts[i] = sendcounts[i];
        self->recvcounts[i] = recvcounts[i];
        self->senddispls[i] = senddispls[i];
        self->recvdispls[i] = recvdispls[i];
    }

    self->nrecv = 0;
    for(i = 0; i < NTask; i ++) {
        self->nrecv += self->recvcounts[i];
    }
    free(sendcounts);
    free(senddispls);
    free(recvcounts);
    free(recvdispls);
    free(ranks);
}

void
fastpm_ghosts_destroy(FastPMGhosts * self)
{
    fastpm_memory_free(self->mem, self->ghostindex);
    free(self->sendcounts);
    free(self->senddispls);
    free(self->recvcounts);
    free(self->recvdispls);
}

static int
 _fastpm_ghosts_probe(FastPMGhosts * self, double x[3], int masterrank, int * ranks)
{
    PM * pm = self->pm;

    int d;
    double low[3];
    double high[3];
    int ilow[3];
    int ihigh[3];
    for(d = 0; d < 3; d ++) {
        low[d] = x[d] - self->margin;
        high[d] = x[d] + self->margin;
    }
    pm_pos_to_rank(pm, low, ilow);
    pm_pos_to_rank(pm, high, ihigh);

    /* probe neighbours */
    int used = 0;
    int j[3];

    for(j[2] = ilow[2]; j[2] <= ihigh[2]; j[2] ++)
    for(j[0] = ilow[0]; j[0] <= ihigh[0]; j[0] ++)
    for(j[1] = ilow[1]; j[1] <= ihigh[1]; j[1] ++) {
        int rank = pm_ipos_to_rank(pm, j);
        if(rank == masterrank)  continue;
        int ptr;
        for(ptr = 0; ptr < used; ptr++) {
            if(rank == ranks[ptr]) break;
        } 
        if(ptr == used) {
            ranks[used++] = rank;
        } 
    }
    return used;
}

void
fastpm_ghosts_fetch(FastPMGhosts * self, FastPMColumn * column, FastPMColumn * recv)
{
    size_t elsize = column->elsize * column->nmemb;

    char * sendbuffer = fastpm_memory_alloc(self->mem, elsize * self->nsend, FASTPM_MEMORY_HEAP);
    char * recvbuffer = fastpm_memory_alloc(self->mem, elsize * self->nrecv, FASTPM_MEMORY_HEAP);
    char * p1 = sendbuffer;
    char * p2 = recvbuffer;

    ptrdiff_t i;
    /* pack the column */

    for(i = 0; i < self->nsend; i ++) {
        fastpm_column_get(column, self->ghostindex[i].i, p1);
        p1 += elsize;
    }

    MPI_Datatype dtype;
    MPI_Type_contiguous(elsize, MPI_BYTE, &dtype);
    MPI_Type_commit(&dtype);
    MPI_Alltoallv(sendbuffer, self->sendcounts, self->senddispls, dtype,
                         recvbuffer, self->recvcounts, self->recvdispls, dtype,
                    self->comm);
    MPI_Type_free(&dtype);

    /* unpack the column */

    fastpm_column_resize(recv, self->nrecv);

    for(i = 0; i < self->nrecv; i ++) {
        fastpm_column_set(recv, i, p2);
        p2 += elsize;
    }

    fastpm_memory_free(self->mem, recvbuffer);
    fastpm_memory_free(self->mem, sendbuffer);
}

void
fastpm_ghosts_reduce(FastPMGhosts * self, FastPMColumn * column, FastPMColumn * local)
{
    size_t elsize = column->elsize * column->nmemb;

    char * sendbuffer = fastpm_memory_alloc(self->mem, elsize * self->nsend, FASTPM_MEMORY_HEAP);
    char * recvbuffer = fastpm_memory_alloc(self->mem, elsize * self->nrecv, FASTPM_MEMORY_HEAP);
    char * p1 = sendbuffer;
    char * p2 = recvbuffer;
    ptrdiff_t i;
    /* pack the column */

    for(i = 0; i < self->nrecv; i ++) {
        fastpm_column_get(local, i, p2);
        p2 += elsize;
    }

    MPI_Datatype dtype;
    MPI_Type_contiguous(elsize, MPI_BYTE, &dtype);
    MPI_Type_commit(&dtype);
    MPI_Alltoallv(recvbuffer, self->recvcounts, self->recvdispls, dtype,
                         sendbuffer, self->sendcounts, self->senddispls, dtype,
                    self->comm);
    MPI_Type_free(&dtype);

    /* unpack the column */

    for(i = 0; i < self->nsend; i ++) {
        double addt[column->nmemb];
        double orig[column->nmemb];

        fastpm_column_to_double(column, p1, addt);
        fastpm_column_get_double(column, self->ghostindex[i].i, orig);
        /* reduction + !*/
        int d;
        for(d = 0; d < column->nmemb; d++) {
            orig[d] += addt[d];
        }
        fastpm_column_set_double(column, self->ghostindex[i].i, orig);
        p1 += elsize;
    }

    fastpm_memory_free(self->mem, recvbuffer);
    fastpm_memory_free(self->mem, sendbuffer);
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    int ThisTask;
    int NTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    FastPMColumn x[1];
    fastpm_column_init_float3(x, 256);

    size_t oldsize = (1 + ThisTask) * 10;
    fastpm_column_resize(x, oldsize);
    int i;

    for(i = 0; i < oldsize; i ++) {
        double pos[3] = {ThisTask * 10 + i, ThisTask * 10 + i, i};
        fastpm_column_set_double(x, i, pos);
    }

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("oldsize = %td\n", oldsize);
    for(i = 0; i < oldsize; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("oldpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }
    //PM pm[1];
    FastPMDomain domain[1];
    //pm_init_simple(pm, 32, 32., comm);
    //fastpm_domain_init(domain, pm, x, comm)

    fastpm_domain_decompose(domain, (FastPMColumn * []) {x, NULL});

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    for(i = 0; i < x->size; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("newpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_destroy(x);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
