#ifndef _FASTPM_STORE_H_
#define _FASTPM_STORE_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS

size_t
fastpm_store_alloc_evenly(FastPMStore * p, size_t np_total, int attributes,
    double alloc_factor, MPI_Comm comm);

void
fastpm_store_alloc(FastPMStore * p, size_t np_upper, int attributes);

void
fastpm_store_init(FastPMStore * p);

void
fastpm_store_destroy(FastPMStore * p);


void
fastpm_store_summary(FastPMStore * p, double dx1[3], double dx2[3], MPI_Comm comm);

void
fastpm_store_set_lagrangian_position(FastPMStore * p, PM * pm, double shift[3], int Nc[3]);

void
fastpm_store_wrap(FastPMStore * p, double BoxSize[3]);

typedef int (fastpm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void
fastpm_store_decompose(FastPMStore * p, fastpm_store_target_func target_func, void * data, MPI_Comm comm);

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

