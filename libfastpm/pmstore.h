#include <stddef.h>
#include <stdint.h>

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

size_t 
pm_store_alloc_evenly(PMStore * p, size_t np_total, int attributes, 
    double alloc_factor, MPI_Comm comm);

void 
pm_store_alloc(PMStore * p, size_t np_upper, int attributes);

void 
pm_store_init(PMStore * p);

void 
pm_store_destroy(PMStore * p);


void 
pm_store_summary(PMStore * p, MPI_Comm comm);

void 
pm_store_set_lagrangian_position(PMStore * p, PM * pm, double shift[3]);

void 
pm_store_wrap(PMStore * p, double BoxSize[3]);

typedef int (pm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void 
pm_store_decompose(PMStore * p, pm_store_target_func target_func, void * data, MPI_Comm comm);

/* Generic IO; unimplemented */
void pm_store_read(PMStore * p, char * datasource);
void pm_store_write(PMStore * p, char * datasource);

void 
pm_store_create_subsample(PMStore * out, PMStore * in, int attributes, int mod, int nc);

