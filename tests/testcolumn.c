#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

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
    void * buffer;

    void
    (*get_double)(FastPMColumn * self, size_t i, double * dest);
    void
    (*set_double)(FastPMColumn * self, size_t i, double * src);

    void
    (* get)(FastPMColumn * self, size_t i, void * dest);
    void
    (* set)(FastPMColumn * self, size_t i, void * src);

    void
    (* destroy)(FastPMColumn * self);

    void * priv;
};


FASTPM_END_DECLS

#endif

#include <mpsort.h>

void
fastpm_column_init(FastPMColumn * self,
        size_t nmemb,
        size_t elsize,
        size_t size)
{
    self->mem = _libfastpm_get_gmem();
    self->elsize = elsize;
    self->nmemb = nmemb;
    self->size = size;
}

static void
fastpm_column_alloc(FastPMColumn * self)
{
    self->buffer = fastpm_memory_alloc(self->mem, self->nmemb * self->elsize * self->size, FASTPM_MEMORY_HEAP);
}
static void
_fastpm_column_destroy(FastPMColumn * self)
{
    if(self->buffer)
        fastpm_memory_free(self->mem, self->buffer);
}

static void *
_fastpm_column_get_dataptr(FastPMColumn * self, size_t i)
{
    return (self->buffer + i * self->elsize * self->nmemb);
}

static void
_fastpm_column_get_double_float64(FastPMColumn * self, size_t i, double * dest)
{
    int d;
    void * ptr = _fastpm_column_get_dataptr(self, i);
    for(d = 0; d < self->nmemb; d ++) {
        dest[d] = ((double*)ptr)[d];
    }
}

static void
_fastpm_column_get_double_float32(FastPMColumn * self, size_t i, double * dest)
{
    int d;
    void * ptr = _fastpm_column_get_dataptr(self, i);
    for(d = 0; d < self->nmemb; d ++) {
        dest[d] = ((float*)ptr)[d];
    }
}

static void
_fastpm_column_get(FastPMColumn * self, size_t i, void * dest)
{
    void * ptr = _fastpm_column_get_dataptr(self, i);
    memcpy(dest, ptr, self->elsize * self->nmemb);
}

static void
_fastpm_column_set_double_float64(FastPMColumn * self, size_t i, double * src)
{
    int d;
    void * ptr = _fastpm_column_get_dataptr(self, i);
    for(d = 0; d < self->nmemb; d ++) {
        ((double*)ptr)[d] = src[d];
    }
}

static void
_fastpm_column_set_double_float32(FastPMColumn * self, size_t i, double * src)
{
    int d;
    void * ptr = _fastpm_column_get_dataptr(self, i);
    for(d = 0; d < self->nmemb; d ++) {
        ((float*)ptr)[d] = src[d];
    }
}

static void
_fastpm_column_set(FastPMColumn * self, size_t i, void * src)
{
    void * ptr = _fastpm_column_get_dataptr(self, i);
    memcpy(ptr, src, self->elsize * self->nmemb);
}

void
fastpm_column_init_double3(FastPMColumn * self, size_t size)
{
    fastpm_column_init(self, 3, sizeof(double), size);
    self->get_double = _fastpm_column_get_double_float64;
    self->set_double = _fastpm_column_set_double_float64;
    self->get = _fastpm_column_get;
    self->set = _fastpm_column_set;
    self->destroy = _fastpm_column_destroy;

    fastpm_column_alloc(self);
}

void
fastpm_column_init_float3(FastPMColumn * self, size_t size)
{
    fastpm_column_init(self, 3, sizeof(float), size);
    self->get_double = _fastpm_column_get_double_float32;
    self->set_double = _fastpm_column_set_double_float32;
    self->get = _fastpm_column_get;
    self->set = _fastpm_column_set;
    self->destroy = _fastpm_column_destroy;
    fastpm_column_alloc(self);
}

void
fastpm_column_init_id(FastPMColumn * self, size_t size)
{
    fastpm_column_init(self, 1, sizeof(uint64_t), size);
    self->get_double = NULL;
    self->set_double = NULL;
    self->get = _fastpm_column_get;
    self->set = _fastpm_column_set;
    self->destroy = _fastpm_column_destroy;
    fastpm_column_alloc(self);
}

void
fastpm_column_destroy(FastPMColumn * self)
{
    self->destroy(self);
}

void
fastpm_column_get_double(FastPMColumn * self, size_t i, double * dest)
{
    self->get_double(self, i, dest);
}

void
fastpm_column_set_double(FastPMColumn * self, size_t i, double * src)
{
    self->set_double(self, i, src);
}

void
fastpm_column_get(FastPMColumn * self, size_t i, double * dest)
{
    self->get(self, i, dest);
}

void
fastpm_column_set(FastPMColumn * self, size_t i, double * src)
{
    self->set(self, i, src);
}

static void
_sort_by_index(const void * ptr, void * radix, void * arg)
{

    int64_t value = *((int64_t *) (ptr));
    value += 9223372036854775808uL;
    *((uint64_t *) radix) = value;
}

void
fastpm_column_parallel_permute(FastPMColumn * self, int64_t index[], MPI_Comm comm)
{
    /* use this to do domain decomposition */

    /* MPSort uses a lot of memory, so we will do one column per time. Thus this api does not do multiple columns */

    char * buffer = fastpm_memory_alloc(self->mem, (sizeof(int64_t) + self->elsize * self->nmemb) * self->size, FASTPM_MEMORY_STACK);
    char * p = buffer;
    /* pack the columns with index */
    size_t packsize = self->elsize * self->nmemb + sizeof(index[0]);
    ptrdiff_t i = 0;
    for(i = 0; i < self->size; i ++) {
        memcpy(p, &index[i], sizeof(index[i]));
        fastpm_column_get(self, i, (void*) (p + sizeof(index[i])));
        p += packsize;
    }

    /* call mpsort */
    mpsort_mpi(buffer, self->size, packsize, _sort_by_index, sizeof(index[0]), p, comm);

    /* unpack the columns */
    p = buffer;
    for(i = 0; i < self->size; i ++) {
        fastpm_column_set(self, i, (void*) (p + sizeof(index[i])));
        p += packsize;
    }
    fastpm_memory_free(self->mem, buffer);
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMColumn x[1];
    int64_t index[128];
    fastpm_column_init_double3(x, 128);
    int i;

    for(i = 0; i < 128; i ++) {
        double pos[3] = {i, i, i};
        fastpm_column_set_double(x, i, pos);
        index[i] = -i;
    }

    fastpm_column_parallel_permute(x, index, comm);

    for(i = 0; i < 128; i +=10) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        fastpm_info("pos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }

    fastpm_column_destroy(x);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
