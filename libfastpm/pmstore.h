#include <stddef.h>
#include <stdint.h>

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
pm_store_summary(PMStore * p, double dx1[3], double dx2[3], MPI_Comm comm);

void 
pm_store_set_lagrangian_position(PMStore * p, PM * pm, double shift[3], int Nc[3]);

void 
pm_store_wrap(PMStore * p, double BoxSize[3]);

typedef int (pm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void 
pm_store_decompose(PMStore * p, pm_store_target_func target_func, void * data, MPI_Comm comm);

/* Generic IO; unimplemented */
void pm_store_read(PMStore * p, char * datasource);
void pm_store_write(PMStore * p, char * datasource);

void 
pm_store_create_subsample(PMStore * out, PMStore * in, int mod, int nc);

