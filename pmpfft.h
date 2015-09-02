#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <pfft.h>

typedef struct {
    void * (*malloc )(size_t);
    void   (*free   )(void *);
    void   (*get_position)(void * pdata, ptrdiff_t index, double pos[3]);
    size_t (*pack)  (void * pdata, ptrdiff_t index, void * packed, int attributes);
    void   (*unpack)(void * pdata, ptrdiff_t index, void * packed, int attributes);
    void   (*reduce)(void * pdata, ptrdiff_t index, void * packed, int method);
} PMIFace;

typedef struct {
    PMIFace iface;

    double (* x)[3];
    float (* v)[3];
    float * acc[3];
    float (* dx1)[3];
    float (* dx2)[3];
    uint64_t * id;
    size_t np;
    size_t np_upper;
} PMStore;

#define PACK_POS   (1 << 0)
#define PACK_VEL   (1 << 1)
#define PACK_DX1   (1 << 2)
#define PACK_DX2   (1 << 3)
#define PACK_ID    (1 << 4)
#define PACK_ACC_X (1 << 5)
#define PACK_ACC_Y (1 << 6)
#define PACK_ACC_Z (1 << 7)
#define HAS(a, b) ((a & b) != 0)

typedef struct {
    ptrdiff_t Nmesh;
    double BoxSize;
    int GhostAttributes;
    int AllAttributes;
} PMInit;

typedef struct {
    ptrdiff_t * edges_int[2];
    double * edges_float[2];
    int * MeshtoCart[2];
} PMGrid;

typedef struct {
    ptrdiff_t start[3];
    ptrdiff_t size[3];
} PMRegion;

typedef struct {
    PMInit init;
    PMIFace iface;
    int NTask;
    int ThisTask;

    pfft_plan r2c;   /* Forward r2c plan */
    pfft_plan c2r;   /* Bacward c2r plan */

    int Nproc[2];
    MPI_Comm Comm2D;

    ptrdiff_t Nmesh[3];
    double    BoxSize[3];

    double    Below[3];
    double    Above[3];

    ptrdiff_t allocsize;

    PMRegion IRegion;
    PMRegion ORegion;
    double * canvas;
    double * workspace;
 
    PMGrid Grid;
    double * MeshtoK[3];

} PM;


typedef struct {
    void * pdata;
    size_t np;
    size_t nghosts;

    /* private members */
    int * Nsend;
    int * Osend;
    int * Nrecv;
    int * Orecv;
    void * send_buffer;
    void * recv_buffer;

    /* iterator status */
    ptrdiff_t ipar; 
    ptrdiff_t ighost;
    int rank;
    int ReductionAttributes;
    size_t elsize;
} PMGhostData;

typedef void (*pm_iter_ghosts_func)(PM * pm, PMGhostData * ppd);


void pm_pfft_init(PM * pm, PMInit * init, PMIFace * iface, MPI_Comm comm);
int pm_pos_to_rank(PM * pm, double pos[3]);

static inline size_t cumsum(int * out, int * in, size_t nitems) {
    size_t total = 0;
    int i;
    for(i = 0; i < nitems; i ++) {
        total += in[i];
        if (out == NULL) continue;
        if(i > 1) 
            out[i] = out[i - 1] + in[i - 1];
        else
            out[i] = 0;
    }
    return total;
}

size_t pm_append_ghosts(PM * pm, size_t np_upper, PMGhostData * ppd);
void pm_reduce_ghosts(PM * pm, PMGhostData * ppd, int attributes);

void pm_paint(PM * pm, void * pdata, ptrdiff_t size);
double pm_readout_one(PM * pm, void * pdata, ptrdiff_t i);
