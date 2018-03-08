#ifndef _FASTPM_STORE_H_
#define _FASTPM_STORE_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS
enum FastPMPackFields {
    PACK_POS   =  1 << 0,
    PACK_VEL   =  1 << 1,
    PACK_DX1   =  1 << 2,
    PACK_DX2   =  1 << 3,
    PACK_ACC   =  1 << 4,
    PACK_ID    =  1 << 5,
    PACK_Q     =  1 << 6,
    PACK_AEMIT     =  1 << 7,
    PACK_POTENTIAL =  1 << 8,
    PACK_TIDAL     =  1 << 9,


    PACK_ACC_X =  1 << 10,
    PACK_ACC_Y =  1 << 11,
    PACK_ACC_Z =  1 << 12,
    PACK_DX1_X   =  1 << 13,
    PACK_DX1_Y   =  1 << 14,
    PACK_DX1_Z   =  1 << 15,
    PACK_DX2_X   =  1 << 16,
    PACK_DX2_Y   =  1 << 17,
    PACK_DX2_Z   =  1 << 18,
    PACK_POS_X =  1 << 19,
    PACK_POS_Y =  1 << 20,
    PACK_POS_Z =  1 << 21,

    PACK_TIDAL_XX =  1 << 22,
    PACK_TIDAL_YY =  1 << 23,
    PACK_TIDAL_ZZ =  1 << 24,
    PACK_TIDAL_XY =  1 << 25,
    PACK_TIDAL_YZ =  1 << 26,
    PACK_TIDAL_ZX =  1 << 27,
};

struct FastPMStore {
    fastpm_posfunc get_position;

    size_t (*pack)  (FastPMStore * p, ptrdiff_t index, void * packed, int attributes);
    void   (*unpack)(FastPMStore * p, ptrdiff_t index, void * packed, int attributes);
    void   (*reduce)(FastPMStore * p, ptrdiff_t index, void * packed, int attributes);
    double (*to_double)(FastPMStore * p, ptrdiff_t index, int attribute);
    void   (*from_double)(FastPMStore * p, ptrdiff_t index, int attribute, double value);

    FastPMMemory * mem;

    int attributes; /* bit flags of allocated attributes */

    double (* x)[3];
    float (* q)[3];
    float (* v)[3];
    float (* acc)[3];
    float (* dx1)[3];
    float (* dx2)[3];
    float (* aemit);
    float (* potential);
    float (* tidal)[6];
    uint64_t * id;
    size_t np;
    size_t np_upper;
    double a_x;
    double a_v;
};

void
fastpm_store_init(FastPMStore * p, size_t np_upper, int attributes, enum FastPMMemoryLocation loc);

size_t
fastpm_store_init_evenly(FastPMStore * p, size_t np_total, int attributes,
    double alloc_factor, MPI_Comm comm);

void
fastpm_store_destroy(FastPMStore * p);

void
fastpm_store_summary(FastPMStore * p, double dx1[3], double dx2[3], MPI_Comm comm);

void
fastpm_store_set_lagrangian_position(FastPMStore * p, PM * pm, double shift[3], ptrdiff_t Nc[3]);

void
fastpm_store_wrap(FastPMStore * p, double BoxSize[3]);

typedef int (fastpm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void
fastpm_store_decompose(FastPMStore * p, fastpm_store_target_func target_func, void * data, MPI_Comm comm);

void fastpm_store_sort_by_id(FastPMStore * p);

size_t
fastpm_store_get_np_total(FastPMStore * p, MPI_Comm comm);

/* Generic IO; unimplemented */
void fastpm_store_read(FastPMStore * p, char * datasource);
void fastpm_store_write(FastPMStore * p, char * datasource);

void
fastpm_store_create_subsample(FastPMStore * out, FastPMStore * in, int mod, int nc);

void
fastpm_store_copy(FastPMStore * in, FastPMStore * out);

void fastpm_store_get_position(FastPMStore * p, ptrdiff_t index, double pos[3]);
void fastpm_store_get_lagrangian_position(FastPMStore * p, ptrdiff_t index, double pos[3]);

FASTPM_END_DECLS

#endif

