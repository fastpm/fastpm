#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <pfft.h>

#define FFT_PRECISION 32
#if FFT_PRECISION == 64
    typedef double real_t;
#elif FFT_PRECISION == 32
    typedef float real_t;
#else
    #error FFT_PRECISION must be 32 or 64
#endif

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

    int attributes; /* bit flags of allocated attributes */

    double (* x)[3];
    float (* q)[3];
    float (* v)[3];
    float (* acc)[3];
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
#define PACK_ACC   (1 << 4)
#define PACK_ID    (1 << 5)
#define PACK_Q     (1 << 6)


#define PACK_ACC_X (1 << 10)
#define PACK_ACC_Y (1 << 11)
#define PACK_ACC_Z (1 << 12)
#define PACK_DX1_X   (1 << 13)
#define PACK_DX1_Y   (1 << 14)
#define PACK_DX1_Z   (1 << 15)
#define PACK_DX2_X   (1 << 16)
#define PACK_DX2_Y   (1 << 17)
#define PACK_DX2_Z   (1 << 18)
#define HAS(a, b) ((a & b) != 0)

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

typedef struct {
    /* in units of real numbers, not bytes. */
    ptrdiff_t start[3];
    ptrdiff_t size[3];
    ptrdiff_t strides[3]; 
    ptrdiff_t total;
} PMRegion;

typedef struct {
    PMInit init;
    PMIFace iface;
    int NTask;
    int ThisTask;

    void * r2c;   /* Forward r2c plan */
    void * c2r;   /* Bacward c2r plan */

    int Nproc[2];
    MPI_Comm Comm2D;

    ptrdiff_t Nmesh[3];
    double    BoxSize[3];

    double    Below[3];
    double    Above[3];

    ptrdiff_t allocsize;
    PMRegion IRegion;
    PMRegion ORegion;
    real_t * canvas;
    real_t * workspace;
 
    PMGrid Grid;
    double * MeshtoK[3];
    double Norm;
    double Volume;
    double CellSize[3];
    double InvCellSize[3];
} PM;


typedef struct {
    PM * pm;
    void * pdata;
    size_t np;
    size_t np_upper;
    size_t nghosts;
    int attributes;
    void   (*get_position)(void * pdata, ptrdiff_t index, double pos[3]);

    /* private members */
    int * Nsend;
    int * Osend;
    int * Nrecv;
    int * Orecv;
    void * send_buffer;
    void * recv_buffer;

    /* iterator status */
    ptrdiff_t ipar; 
    int * ighost_to_ipar;
    int rank;
    ptrdiff_t * reason; /* relative offset causing the ghost */
    int ReductionAttributes;
    size_t elsize;
} PMGhostData;

typedef void (*pm_iter_ghosts_func)(PM * pm, PMGhostData * ppd);

typedef struct {
    float k_finite; /* k, finite */
    float k; /* k */
    float kk_finite; /* k ** 2, on a mesh */
    float kk;  /* k ** 2 */
    float cic;  /* 1 - 2 / 3 sin^2 ( 0.5 k L / N)*/
    float extra;  /* any temporary variable that can be useful. */
} PMKFactors;


#ifdef __cplusplus
extern "C" {
#endif
void pm_init(PM * pm, PMInit * init, PMIFace * iface, MPI_Comm comm);

void 
pm_init_simple(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void pm_destroy(PM * pm);
void pm_start(PM * pm);
void pm_stop(PM * pm);

int pm_pos_to_rank(PM * pm, double pos[3]);
int pm_ipos_to_rank(PM * pm, int i[3]);

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

void pm_unravel_o_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]);
void pm_inc_o_index(PM * pm, ptrdiff_t i[3]);
void pm_unravel_i_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]);
void pm_inc_i_index(PM * pm, ptrdiff_t i[3]);

void pm_append_ghosts(PMGhostData * pgd);
void pm_reduce_ghosts(PMGhostData * pgd, int attributes);
void pm_destroy_ghosts(PMGhostData * pgd);

void pm_paint(PM * pm, void * pdata, ptrdiff_t size);
double pm_readout_pos(PM * pm, double pos[3]);
void pm_paint_pos(PM * pm, double pos[3], double weight);
double pm_readout_one(PM * pm, PMStore * p, ptrdiff_t i);
#ifdef __cplusplus
}
#endif

typedef int (pm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void pm_store_read(PMStore * p, char * datasource);
void pm_store_write(PMStore * p, char * datasource);
void pm_store_destroy(PMStore * p);

void pm_store_init(PMStore * p);
void pm_store_alloc(PMStore * p, size_t np_upper, int attributes);

size_t 
pm_store_alloc_evenly(PMStore * p, size_t np_total, int attributes, double alloc_factor, MPI_Comm comm);
void pm_store_decompose(PMStore * p, pm_store_target_func target_func, void * data, MPI_Comm comm);
void pm_store_wrap(PMStore * p, double BoxSize[3]);

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);

void 
pm_create_k_factors(PM * pm, PMKFactors * fac[3]);
void 
pm_destroy_k_factors(PM * pm, PMKFactors * fac[3]);

int omp_get_num_threads();
int omp_get_thread_num();

static inline
void
pm_prepare_omp_loop(PM * pm, ptrdiff_t * start, ptrdiff_t * end, ptrdiff_t i[3]) 
{ 
    /* static schedule the openmp loops. start, end is in units of 'real' numbers.
     *
     * i is in units of complex numbers.
     *
     * We call pm_unravel_o_index to set the initial i[] for each threads,
     * then rely on pm_inc_o_index to increment i, because the former is 
     * much slower than pm_inc_o_index and would eliminate threading advantage.
     *
     * */
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();

    *start = ith * pm->ORegion.total / nth * 2;
    *end = (ith + 1) * pm->ORegion.total / nth * 2;

    /* do not unravel if we are not looping at all. 
     * This fixes a FPE when
     * the rank has ORegion.total == 0 
     * -- with PFFT the last transposed dimension
     * on some ranks will be 0 */
    if(*end > *start) 
        pm_unravel_o_index(pm, *start / 2, i);

#if 0
        msg_aprintf(info, "ith %d nth %d start %td end %td pm->ORegion.strides = %td %td %td\n", ith, nth,
            *start, *end,
            pm->ORegion.strides[0],
            pm->ORegion.strides[1],
            pm->ORegion.strides[2]
            );
#endif

}

#ifdef __cplusplus
extern "C" {
#endif
void 
pm_store_set_lagrangian_position(PMStore * p, PM * pm, double shift[3]);

void 
pm_calculate_forces(PMStore * p, PM * pm, double density_factor);

void 
pm_r2c(PM * pm);

void 
pm_c2r(PM * pm);

void 
pm_get_times(int istep,
    double time_step[],
    int nstep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1);
#ifdef __cplusplus
}
#endif
