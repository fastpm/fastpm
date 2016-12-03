FASTPM_BEGIN_DECLS

typedef struct FastPMDomain {
    FastPMMemory * mem;
    FastPMMesh * mesh;
    ptrdiff_t * masterindex;
    MPI_Comm comm;
    size_t oldsize;
    size_t newsize;
} FastPMDomain;


void 
fastpm_domain_init(FastPMDomain * self, FastPMMesh * mesh, FastPMColumn * pos, MPI_Comm comm);

void
fastpm_domain_decompose(FastPMDomain * self, FastPMColumn * columns[]);

void
fastpm_domain_destroy(FastPMDomain * self);

typedef struct FastPMGhosts {
    FastPMMemory * mem;
    FastPMMesh * mesh;
    struct FastPMGhostPair * ghostindex;
    MPI_Comm comm;

    size_t nsend;
    size_t nrecv;
    int * sendcounts;
    int * senddispls;
    int * recvcounts;
    int * recvdispls;
} FastPMGhosts;

void
fastpm_ghosts_init(FastPMGhosts * self, FastPMMesh * mesh, FastPMColumn * pos, FastPMColumn * margin, MPI_Comm comm);

void
fastpm_ghosts_destroy(FastPMGhosts * self);

FastPMColumn *
fastpm_ghosts_fetch(FastPMGhosts * self, FastPMColumn * column);

void
fastpm_ghosts_reduce(FastPMGhosts * self, FastPMColumn * column, FastPMColumn * local);

FASTPM_END_DECLS
