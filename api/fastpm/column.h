#ifndef _FASTPM_COLUMN_H_
#define _FASTPM_COLUMN_H_

#include <stddef.h>
#include <stdint.h>

FASTPM_BEGIN_DECLS

typedef struct FastPMColumn FastPMColumn;

struct FastPMColumn {
    FastPMMemory * mem;
    size_t nmemb;
    size_t elsize;
    size_t size;
    size_t maxsize;
    void * buffer;

    void
    (*to_double)(FastPMColumn * self, void * ptr, double * dest);
    void
    (*from_double)(FastPMColumn * self, void * ptr, double * src);

    void
    (* get)(FastPMColumn * self, ptrdiff_t i, void * dest);
    void
    (* set)(FastPMColumn * self, ptrdiff_t i, void * src);

    void
    (* destroy)(FastPMColumn * self);

    void
    (* resize)(FastPMColumn * self, size_t newsize);

    void * priv;

    double timestamp; /* a */
};

void
fastpm_column_init(FastPMColumn * self,
        size_t nmemb,
        size_t elsize,
        size_t maxsize);

void
fastpm_column_resize(FastPMColumn * self, size_t size);

void
fastpm_column_init_double3(FastPMColumn * self, size_t maxsize);

void
fastpm_column_init_float3(FastPMColumn * self, size_t maxsize);

void
fastpm_column_init_uint64(FastPMColumn * self, size_t maxsize);

void
fastpm_column_destroy(FastPMColumn * self);

void
fastpm_column_get_double(FastPMColumn * self, ptrdiff_t i, double * dest);

void
fastpm_column_set_double(FastPMColumn * self, ptrdiff_t i, double * src);

void
fastpm_column_to_double(FastPMColumn * self, void * ptr, double * dest);

void
fastpm_column_from_double(FastPMColumn * self, void *ptr, double * src);

void
fastpm_column_get(FastPMColumn * self, ptrdiff_t i, void * dest);

void
fastpm_column_set(FastPMColumn * self, ptrdiff_t i, void * src);

void
fastpm_column_parallel_permute(FastPMColumn * self, int64_t index[], size_t newsize, MPI_Comm comm);

FASTPM_END_DECLS

#endif
