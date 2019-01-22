#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>

#ifdef ENABLE_VALGRIND
#include </usr/include/valgrind/memcheck.h>
#else
#define VALGRIND_MAKE_MEM_DEFINED(a, b) do {} while(0);
#define VALGRIND_CHECK_MEM_IS_DEFINED(a, b) do {} while(0);
#endif


#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#ifndef LIKELY
#define LIKELY(x) __builtin_expect(!!(x), 1)
#endif

#ifndef UNLIKELY
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

typedef struct {
    ptrdiff_t Nmesh;
    double BoxSize;
    int NprocY;
    int transposed;
    int use_fftw;
} PMInit;

typedef struct {
    ptrdiff_t * edges_int[2];
    double * edges_float[2];
    int * MeshtoCart[2];
} PMGrid;

struct PM {
    PMInit init;

    int NTask;
    int ThisTask;

    void * r2c;   /* Forward r2c plan */
    void * c2r;   /* Bacward c2r plan */

    int Nproc[2];
    MPI_Comm Comm2D;

    ptrdiff_t Nmesh[3];
    double    BoxSize[3];

    ptrdiff_t allocsize;
    PMRegion IRegion;
    PMRegion ORegion;

    PMGrid Grid;
    double * MeshtoK[3];
    double Norm;
    double Volume;
    double CellSize[3];
    double InvCellSize[3];

    FastPMMemory * mem;
};

void
pm_module_init();

void
pm_module_cleanup();

/* Initializing a PM object. */
void 
pm_init(PM * pm, PMInit * init, MPI_Comm comm);

void 
pm_init_simple(PM * pm, int Ngrid, double BoxSize, MPI_Comm comm);

void pm_destroy(PM * pm);

int pm_pos_to_rank(PM * pm, double pos[3]);
int pm_ipos_to_rank(PM * pm, int i[3]);

/* reset 'x' and 'q' of every particle to the lagrangian position. This function shall
 * not belong here.*/

/* This function guarentees good performance for sparse all to all. 
 * Used e.g. in domain decomposition */
int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);

static inline size_t cumsum(int * out, int * in, size_t nitems) {
    size_t total = 0;
    int i;
    for(i = 0; i < nitems; i ++) {
        total += in[i];
        if (out == NULL) continue;
        if(i >= 1) 
            out[i] = out[i - 1] + in[i - 1];
        else
            out[i] = 0;
    }
    return total;
}

