#ifndef _FASTPM_STORE_H_
#define _FASTPM_STORE_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS

typedef uint8_t FastPMParticleMaskType;

/*
 * A solver has multiple particle species;
 * the species names are defined here.
 * */
enum FastPMSpecies {
    FASTPM_SPECIES_BARYON = 0,
    FASTPM_SPECIES_CDM = 1,
    FASTPM_SPECIES_NCDM = 2,
};

const char *
fastpm_species_get_name(enum FastPMSpecies species);

double
fastpm_store_get_mass(FastPMStore * p, ptrdiff_t index);

/*
 * enum constants for naming the columns in the source code.
 * Keep the value ordering in agreement with the ordering in FastPMStore struct;
 * if mismatched, the code will issue an error at start up.
 * */
typedef enum FastPMColumnTags {
    COLUMN_MASK    =  1L << 0,
    COLUMN_POS   =  1L << 1,
    COLUMN_Q     =  1L << 2,
    COLUMN_VEL   =  1L << 3,
    COLUMN_DX1   =  1L << 4,
    COLUMN_DX2   =  1L << 5,
    COLUMN_DV1   =  1L << 6,
    COLUMN_ACC   =  1L << 7,
    COLUMN_ID    =  1L << 8,
    COLUMN_AEMIT     =  1L << 9,
    COLUMN_DENSITY =  1L << 10,
    COLUMN_POTENTIAL =  1L << 11,
    COLUMN_TIDAL     =  1L << 12,
    COLUMN_PGDC   =  1L << 13,

    /* for fof */
    COLUMN_MINID =  1L << 14,
    COLUMN_TASK =  1L << 15,
    COLUMN_LENGTH =  1L << 16,
    COLUMN_RDISP =  1L << 17,
    COLUMN_VDISP =  1L << 18,
    COLUMN_RVDISP =  1L << 19,
    
    COLUMN_MASS = 1L << 20,

} FastPMColumnTags;

struct FastPMStore {
    FastPMMemory * mem;
    char name[32];
    FastPMColumnTags attributes; /* bit flags of allocated attributes */

    void * _base; /* base pointer of all memory buffers */

    size_t np;
    size_t np_upper;

    /* The ordering of the column_info array is the same as the columns array */
    struct FastPMColumnInfo {
        void   (*pack)   (FastPMStore * p, ptrdiff_t index, int ci, void * packed);
        void   (*unpack) (FastPMStore * p, ptrdiff_t index, int ci, void * packed);
        double (*to_double) (FastPMStore * p, ptrdiff_t index, int ci, int memb);
        void   (*from_double) (FastPMStore * p, ptrdiff_t index, int ci, int memb, const double value);

        char name[32];
        char dtype[8];
        size_t elsize;
        size_t membsize;
        size_t nmemb;

        /* This is exactly 1 << item index. Ensured by DEFINE_COLUMN and values in FastPMColumnTags. */
        FastPMColumnTags attribute;
    } _column_info[32];

    /* meta-data: remember to modify io.c after adding an item. */
    struct {
        double a_x;
        double a_v;
        double M0;    /* base mass in 10^10 M_sun / h; particle mass is M0 + mass[i] */

        double _q_shift[3];
        double _q_scale[3];
        ptrdiff_t _q_strides[3];
        ptrdiff_t _q_size;
    } meta;

    union {
        char * columns[32];
        struct {
            FastPMParticleMaskType * mask;
            double (* x)[3];
            float (* q)[3];
            float (* v)[3];  
            float (* dx1)[3];
            float (* dx2)[3];
            float (* dv1)[3];
            float (* acc)[3];
            uint64_t * id;
            float (* aemit);
            float (* rho);
            float (* potential);
            float (* tidal)[6];
            float (* pgdc)[3];

            /* for fof */
            uint64_t * minid;
            int32_t  * task;
            int32_t * length;
            float (* rdisp)[6]; /* zero lag, first lag, second lag */
            float (* vdisp)[6];
            float (* rvdisp)[9];

            /* multiple species support */
            float (* mass);   /* extra mass in addition to meta.M0; see fastpm_store_get_mass */
        };
    };
};

/* convert a column name literal (e.g. x, id) to the column index in the column_info array */
#define FASTPM_STORE_COLUMN_INDEX(column) (((char*) &(((FastPMStore *) NULL)->column) - (char*) &(((FastPMStore *)NULL)->columns[0])) \
                        / sizeof(((FastPMStore *) NULL)->columns[0]))
#define FASTPM_STORE_COLUMN_INFO(p, column) ((p)->_column_info[FASTPM_STORE_COLUMN_INDEX(column)])

/* */
typedef struct {
    int elsize;
    int Ncolumns;
    FastPMColumnTags attributes;

    /* private : */
    int _ci[32];
    int _offsets[32];
    struct FastPMColumnInfo _column_info[32];
} FastPMPackingPlan;

typedef struct {
    FastPMColumnTags attribute;
    int memb;   // for a vector field, this number represents the component
} FastPMFieldDescr;

const static FastPMFieldDescr FASTPM_FIELD_DESCR_NONE = {0, 0};
void
fastpm_store_init_details(FastPMStore * p, const char * name, size_t np_upper, FastPMColumnTags attributes, enum FastPMMemoryLocation loc, const char * file, const int line);

#define fastpm_store_init(p, name, np_upper, attributes, loc) fastpm_store_init_details(p, name, np_upper, attributes, loc, __FILE__, __LINE__)

void
fastpm_store_set_name(FastPMStore * store, const char * name);

#define fastpm_store_init_evenly(p, name, np_total, attributes, alloc_factor, comm) fastpm_store_init_evenly_details(p, name, np_total, attributes, alloc_factor, comm, __FILE__, __LINE__)
size_t
fastpm_store_init_evenly_details(FastPMStore * p, const char * name, size_t np_total, FastPMColumnTags attributes,
    double alloc_factor, MPI_Comm comm, const char * file, const int line);

void
fastpm_store_fill(FastPMStore * p, PM * pm, double * shift, ptrdiff_t * Nc);

int
fastpm_store_has_q(FastPMStore *q);

void
fastpm_store_get_q_from_id(FastPMStore * p, uint64_t id, double q[3]);

void
fastpm_store_get_iq_from_id(FastPMStore * p, uint64_t id, ptrdiff_t pabs[3]);

size_t fastpm_store_pack   (FastPMStore * p, ptrdiff_t index, void * packed, FastPMColumnTags attributes);
void   fastpm_store_unpack (FastPMStore * p, ptrdiff_t index, void * packed, FastPMColumnTags attributes);

int fastpm_store_find_column_id(FastPMStore *p, FastPMColumnTags attribute);

void
fastpm_packing_plan_init(FastPMPackingPlan * plan, FastPMStore * p, FastPMColumnTags attributes);

void
fastpm_packing_plan_pack(FastPMPackingPlan * plan,
            FastPMStore * p, ptrdiff_t i, void * packed);

void
fastpm_packing_plan_unpack(FastPMPackingPlan * plan,
            FastPMStore * p, ptrdiff_t i, void * packed);

void
fastpm_packing_plan_unpack_ci(FastPMPackingPlan * plan, int ci,
            FastPMStore * p, ptrdiff_t i, void * packed);

void
fastpm_store_destroy(FastPMStore * p);

void
fastpm_store_summary(FastPMStore * p,
        FastPMColumnTags attribute,
        MPI_Comm comm,
        const char * fmt,
        ...);

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
FastPMReduceAddFloat(FastPMStore * src, ptrdiff_t isrc, FastPMStore * dest, ptrdiff_t idest, int ci, void * userdata);

void
fastpm_store_sort(FastPMStore * p, int (*cmpfunc)(const int i1, const int i2, FastPMStore * p));

size_t
fastpm_store_get_np_total(FastPMStore * p, MPI_Comm comm);

size_t
fastpm_store_get_mask_sum(FastPMStore * p, MPI_Comm comm);

void
fastpm_store_fill_subsample_mask(FastPMStore * p,
        double fraction,
        FastPMParticleMaskType * mask,
        MPI_Comm comm);

void
fastpm_store_fill_subsample_mask_every_dim(FastPMStore * p,
                                              int every,
                                              FastPMParticleMaskType * mask);

size_t
fastpm_store_subsample(FastPMStore * in, FastPMParticleMaskType * mask, FastPMStore * out);

void
fastpm_store_histogram_aemit_sorted(FastPMStore * store,
        int64_t * hist,
        double * edges,
        size_t nbins,
        MPI_Comm comm);

void
fastpm_store_copy(FastPMStore * in, FastPMStore * out);

void
fastpm_store_steal(FastPMStore * in, FastPMStore * out, FastPMColumnTags attributes);

void
fastpm_store_take(FastPMStore * in, ptrdiff_t i, FastPMStore * out, ptrdiff_t j);

void
fastpm_store_extend(FastPMStore * p, FastPMStore * extra);

void fastpm_store_get_position(FastPMStore * p, ptrdiff_t index, double pos[3]);
void fastpm_store_get_lagrangian_position(FastPMStore * p, ptrdiff_t index, double pos[3]);

int
FastPMTargetPM (FastPMStore * p, ptrdiff_t i, PM * pm);

FASTPM_END_DECLS

#endif

