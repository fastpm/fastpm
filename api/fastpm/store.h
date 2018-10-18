#ifndef _FASTPM_STORE_H_
#define _FASTPM_STORE_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS
enum FastPMPackFields {
    PACK_MASK    =  1L << 0,
    PACK_POS   =  1L << 1,
    PACK_Q     =  1L << 2,
    PACK_VEL   =  1L << 3,
    PACK_DX1   =  1L << 4,
    PACK_DX2   =  1L << 5,
    PACK_ACC   =  1L << 6,
    PACK_ID    =  1L << 7,
    PACK_AEMIT     =  1L << 8,
    PACK_DENSITY =  1L << 9,
    PACK_POTENTIAL =  1L << 10,
    PACK_TIDAL     =  1L << 11,

    /* for fof */
    PACK_MINID =  1L << 12,
    PACK_TASK =  1L << 13,
    PACK_LENGTH =  1L << 14,
    PACK_RDISP =  1L << 15,
    PACK_VDISP =  1L << 16,
    PACK_RVDISP =  1L << 17,
};

struct FastPMStore {
    FastPMMemory * mem;

    enum FastPMPackFields attributes; /* bit flags of allocated attributes */

    void * base; /* base pointer of all memory buffers */

    struct {
        void   (*pack)   (FastPMStore * p, ptrdiff_t index, int ci, void * packed);
        void   (*unpack) (FastPMStore * p, ptrdiff_t index, int ci, void * packed);
        void   (*pack_member)   (FastPMStore * p, ptrdiff_t index, int ci, int memb, void * packed);
        void   (*reduce_member) (FastPMStore * p, ptrdiff_t index, int ci, int memb, void * packed);
        double (*to_double) (FastPMStore * p, ptrdiff_t index, int ci, int memb);
        void   (*from_double) (FastPMStore * p, ptrdiff_t index, int ci, int memb, const double value);

        char dtype[8];
        size_t elsize;
        size_t membsize;
        size_t nmemb;
        enum FastPMPackFields attribute;
    } column_info[32];

    union {
        char * columns[32];
        struct {
            uint8_t * mask;
            double (* x)[3];
            float (* q)[3];
            float (* v)[3];
            float (* dx1)[3];
            float (* dx2)[3];
            float (* acc)[3];
            uint64_t * id;
            float (* aemit);
            float (* rho);
            float (* potential);
            float (* tidal)[6];

            /* for fof */
            uint64_t * minid;
            int32_t  * task;
            int32_t * length;
            float (* rdisp)[6]; /* zero lag, first lag, second lag */
            float (* vdisp)[6];
            float (* rvdisp)[9];
        };
    };
    size_t np;
    size_t np_upper;
    double a_x;
    double a_v;

    double q_shift[3];
    double q_scale[3];
    ptrdiff_t q_strides[3];
};

typedef struct {
    enum FastPMPackFields attribute;
    int memb;
} FastPMFieldDescr;

const static FastPMFieldDescr FASTPM_FIELD_DESCR_NONE = {0, 0};
void
fastpm_store_init_details(FastPMStore * p, size_t np_upper, enum FastPMPackFields attributes, enum FastPMMemoryLocation loc, const char * file, const int line);

#define fastpm_store_init(p, np_upper, attributes, loc) fastpm_store_init_details(p, np_upper, attributes, loc, __FILE__, __LINE__)

size_t
fastpm_store_init_evenly(FastPMStore * p, size_t np_total, enum FastPMPackFields attributes,
    double alloc_factor, MPI_Comm comm);

void
fastpm_store_fill(FastPMStore * p, PM * pm, double * shift, ptrdiff_t * Nc);

void 
fastpm_store_get_q_from_id(FastPMStore * p, uint64_t id, double q[3]);

size_t fastpm_store_pack   (FastPMStore * p, ptrdiff_t index, void * packed, enum FastPMPackFields attributes);
void   fastpm_store_unpack (FastPMStore * p, ptrdiff_t index, void * packed, enum FastPMPackFields attributes);

int fastpm_store_find_column_id(FastPMStore *p, enum FastPMPackFields attribute);

void
fastpm_store_destroy(FastPMStore * p);

void
fastpm_store_summary(FastPMStore * p, double dx1[3], double dx2[3], MPI_Comm comm);

void
fastpm_store_wrap(FastPMStore * p, double BoxSize[3]);

typedef int (*fastpm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

int
fastpm_store_decompose(FastPMStore * p, fastpm_store_target_func target_func, void * data, MPI_Comm comm);

void
fastpm_store_permute(FastPMStore * p, int * ind);

int
FastPMLocalSortByID(const int i1,
                    const int i2,
                    FastPMStore * p);

void
fastpm_store_sort(FastPMStore * p, int (*cmpfunc)(const int i1, const int i2, FastPMStore * p));

size_t
fastpm_store_get_np_total(FastPMStore * p, MPI_Comm comm);

size_t
fastpm_store_get_mask_sum(FastPMStore * p, MPI_Comm comm);

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
fastpm_store_histogram_aemit_sorted(FastPMStore * store,
        int64_t * hist,
        double * edges,
        size_t nbins,
        MPI_Comm comm);

void
fastpm_store_copy(FastPMStore * in, FastPMStore * out);

void
fastpm_store_take(FastPMStore * in, ptrdiff_t i, FastPMStore * out, ptrdiff_t j);

void
fastpm_store_append(FastPMStore * in, FastPMStore * out);

void fastpm_store_get_position(FastPMStore * p, ptrdiff_t index, double pos[3]);
void fastpm_store_get_lagrangian_position(FastPMStore * p, ptrdiff_t index, double pos[3]);

int
FastPMTargetPM (FastPMStore * p, ptrdiff_t i, PM * pm);

FASTPM_END_DECLS

#endif

