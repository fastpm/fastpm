#ifndef _FASTPM_STORE_H_
#define _FASTPM_STORE_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS
enum FastPMPackFields {
    PACK_POS   =  1L << 0,
    PACK_VEL   =  1L << 1,
    PACK_DX1   =  1L << 2,
    PACK_DX2   =  1L << 3,
    PACK_ACC   =  1L << 4,
    PACK_ID    =  1L << 5,
    PACK_Q     =  1L << 6,
    PACK_AEMIT     =  1L << 7,
    PACK_DENSITY =  1L << 8,
    PACK_POTENTIAL =  1L << 9,
    PACK_TIDAL     =  1L << 10,
    PACK_FOF =  1L << 11,
    PACK_LENGTH =  1L << 12,
    PACK_MASK    =  1L << 13,


    PACK_ACC_X =  1L << 20,
    PACK_ACC_Y =  1L << 21,
    PACK_ACC_Z =  1L << 22,
    PACK_DX1_X   =  1L << 23,
    PACK_DX1_Y   =  1L << 24,
    PACK_DX1_Z   =  1L << 25,
    PACK_DX2_X   =  1L << 26,
    PACK_DX2_Y   =  1L << 27,
    PACK_DX2_Z   =  1L << 28,
    PACK_POS_X =  1L << 29,
    PACK_POS_Y =  1L << 30,
    PACK_POS_Z =  1L << 31,

    PACK_TIDAL_XX =  1L << 32,
    PACK_TIDAL_YY =  1L << 33,
    PACK_TIDAL_ZZ =  1L << 34,
    PACK_TIDAL_XY =  1L << 35,
    PACK_TIDAL_YZ =  1L << 36,
    PACK_TIDAL_ZX =  1L << 37,
};

struct FastPMStore {
    fastpm_posfunc get_position;

    size_t (*pack)  (FastPMStore * p, ptrdiff_t index, void * packed, enum FastPMPackFields attributes);
    void   (*unpack)(FastPMStore * p, ptrdiff_t index, void * packed, enum FastPMPackFields attributes);
    void   (*reduce)(FastPMStore * p, ptrdiff_t index, void * packed, enum FastPMPackFields attributes);
    double (*to_double)(FastPMStore * p, ptrdiff_t index, enum FastPMPackFields attribute);
    void   (*from_double)(FastPMStore * p, ptrdiff_t index, enum FastPMPackFields attribute, double value);

    FastPMMemory * mem;

    enum FastPMPackFields attributes; /* bit flags of allocated attributes */

    double (* x)[3];
    float (* q)[3];
    float (* v)[3];
    float (* acc)[3];
    float (* dx1)[3];
    float (* dx2)[3];
    float (* aemit);
    float (* rho);
    float (* potential);
    float (* tidal)[6];
    uint64_t * id;
    uint8_t * mask;

    /* for fof */
    struct FastPMFOFData {
        uint64_t minid;
        int64_t task; /* to fill up the 8 bytes alignment */
     } * fof;
    int32_t * length;

    size_t np;
    size_t np_upper;
    double a_x;
    double a_v;

    double q_shift[3];
    double q_scale[3];
    ptrdiff_t q_strides[3];
};

void
fastpm_store_init(FastPMStore * p, size_t np_upper, enum FastPMPackFields attributes, enum FastPMMemoryLocation loc);

size_t
fastpm_store_init_evenly(FastPMStore * p, size_t np_total, enum FastPMPackFields attributes,
    double alloc_factor, MPI_Comm comm);

void
fastpm_store_fill(FastPMStore * p, PM * pm, double * shift, ptrdiff_t * Nc);

void 
fastpm_store_get_q_from_id(FastPMStore * p, uint64_t id, double q[3]);

void
fastpm_store_destroy(FastPMStore * p);

void
fastpm_store_summary(FastPMStore * p, double dx1[3], double dx2[3], MPI_Comm comm);

void
fastpm_store_wrap(FastPMStore * p, double BoxSize[3]);

typedef int (*fastpm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void
fastpm_store_decompose(FastPMStore * p, fastpm_store_target_func target_func, void * data, MPI_Comm comm);

void fastpm_store_sort_by_id(FastPMStore * p);

size_t
fastpm_store_get_np_total(FastPMStore * p, MPI_Comm comm);

/* Generic IO; unimplemented */
void fastpm_store_read(FastPMStore * p, char * datasource);
void fastpm_store_write(FastPMStore * p, char * datasource);

void
fastpm_store_fill_subsample_mask(FastPMStore * p,
        double fraction,
        uint8_t * mask,
        MPI_Comm comm);

void
fastpm_store_subsample(FastPMStore * in, uint8_t * mask, FastPMStore * out);

void
fastpm_store_histogram_aemit(FastPMStore * store,
        ptrdiff_t * hist,
        double * edges,
        size_t nbins,
        MPI_Comm comm);

void
fastpm_store_copy(FastPMStore * in, FastPMStore * out);

void
fastpm_store_append(FastPMStore * in, FastPMStore * out);

void fastpm_store_get_position(FastPMStore * p, ptrdiff_t index, double pos[3]);
void fastpm_store_get_lagrangian_position(FastPMStore * p, ptrdiff_t index, double pos[3]);

int
FastPMTargetPM (FastPMStore * p, ptrdiff_t i, PM * pm);

FASTPM_END_DECLS

#endif

