#include <mpi.h>
#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#include <mpsort.h>

void
fastpm_column_init(FastPMColumn * self,
        size_t nmemb,
        size_t elsize,
        size_t maxsize)
{
    self->mem = _libfastpm_get_gmem();
    self->elsize = elsize;
    self->nmemb = nmemb;
    self->size = 0;
    self->maxsize = maxsize;
}

void
fastpm_column_resize(FastPMColumn * self, size_t size)
{
    if(size >= self->maxsize) {
        fastpm_raise(-1, "requested size exceed the max (%td > %td)\n", size, self->maxsize);
    }
    self->size = size;
}

static void
fastpm_column_alloc(FastPMColumn * self)
{
    self->buffer = fastpm_memory_alloc(self->mem, self->nmemb * self->elsize * self->maxsize, FASTPM_MEMORY_HEAP);
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
_fastpm_column_to_double_float64(FastPMColumn * self, void * ptr, double * dest)
{
    int d;
    for(d = 0; d < self->nmemb; d ++) {
        dest[d] = ((double*)ptr)[d];
    }
}

static void
_fastpm_column_to_double_float32(FastPMColumn * self, void * ptr, double * dest)
{
    int d;
    for(d = 0; d < self->nmemb; d ++) {
        dest[d] = ((float*)ptr)[d];
    }
}

static void
_fastpm_column_get(FastPMColumn * self, ptrdiff_t i, void * dest)
{
    void * ptr = _fastpm_column_get_dataptr(self, i);
    memcpy(dest, ptr, self->elsize * self->nmemb);
}

static void
_fastpm_column_from_double_float64(FastPMColumn * self, void * ptr, double * src)
{
    int d;
    for(d = 0; d < self->nmemb; d ++) {
        ((double*)ptr)[d] = src[d];
    }
}

static void
_fastpm_column_from_double_float32(FastPMColumn * self, void * ptr, double * src)
{
    int d;
    for(d = 0; d < self->nmemb; d ++) {
        ((float*)ptr)[d] = src[d];
    }
}

static void
_fastpm_column_set(FastPMColumn * self, ptrdiff_t i, void * src)
{
    void * ptr = _fastpm_column_get_dataptr(self, i);
    memcpy(ptr, src, self->elsize * self->nmemb);
}

void
fastpm_column_init_double3(FastPMColumn * self, size_t maxsize)
{
    fastpm_column_init(self, 3, sizeof(double), maxsize);
    self->to_double = _fastpm_column_to_double_float64;
    self->from_double = _fastpm_column_from_double_float64;
    self->get = _fastpm_column_get;
    self->set = _fastpm_column_set;
    self->destroy = _fastpm_column_destroy;

    fastpm_column_alloc(self);
}

void
fastpm_column_init_float3(FastPMColumn * self, size_t maxsize)
{
    fastpm_column_init(self, 3, sizeof(float), maxsize);
    self->to_double = _fastpm_column_to_double_float32;
    self->from_double = _fastpm_column_from_double_float32;
    self->get = _fastpm_column_get;
    self->set = _fastpm_column_set;
    self->destroy = _fastpm_column_destroy;
    fastpm_column_alloc(self);
}

void
fastpm_column_init_uint64(FastPMColumn * self, size_t maxsize)
{
    fastpm_column_init(self, 1, sizeof(uint64_t), maxsize);
    self->to_double = NULL;
    self->from_double = NULL;
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
fastpm_column_to_double(FastPMColumn * self, void * ptr, double * dest)
{
    self->to_double(self, ptr, dest);
}

void
fastpm_column_from_double(FastPMColumn * self, void * ptr, double * src)
{
    self->from_double(self, ptr, src);
}

void
fastpm_column_get_double(FastPMColumn * self, ptrdiff_t i, double * dest)
{
    char ptr[self->elsize * self->nmemb];
    fastpm_column_get(self, i, ptr);
    self->to_double(self, ptr, dest);
}

void
fastpm_column_set_double(FastPMColumn * self, ptrdiff_t i, double * src)
{
    char ptr[self->elsize * self->nmemb];
    self->from_double(self, ptr, src);
    fastpm_column_set(self, i, ptr);
}

void
fastpm_column_get(FastPMColumn * self, ptrdiff_t i, void * dest)
{
    self->get(self, i, dest);
}

void
fastpm_column_set(FastPMColumn * self, ptrdiff_t i, void * src)
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
fastpm_column_parallel_permute(FastPMColumn * self, int64_t index[], size_t newsize, MPI_Comm comm)
{
    /* MPSort uses a lot of memory, so we will do one column per time. Thus this api does not do multiple columns */

    /* use this to do domain decomposition */

    /* sanity checks */
    long long sizesum[2] = {newsize, self->size};
    MPI_Allreduce(MPI_IN_PLACE, sizesum, 2, MPI_LONG_LONG, MPI_SUM, comm);
    if (sizesum[0] != sizesum[1]) {
        fastpm_raise(-1, "Number of items is not conserved %ld(new) != %ld(old)\n", sizesum[0], sizesum[1]);
    }

    size_t oldsize = self->size;
    fastpm_column_resize(self, newsize);

    /* pack the columns with index */
    size_t packsize = self->elsize * self->nmemb + sizeof(index[0]);
    size_t buffersize = (oldsize > newsize) ? oldsize : newsize;
    char * buffer = fastpm_memory_alloc(self->mem, (sizeof(int64_t) + self->elsize * self->nmemb) * buffersize, FASTPM_MEMORY_STACK);
    char * p = buffer;

    ptrdiff_t i = 0;
    for(i = 0; i < oldsize; i ++) {
        memcpy(p, &index[i], sizeof(index[i]));
        fastpm_column_get(self, i, (void*) (p + sizeof(index[i])));
        p += packsize;
    }

    /* call mpsort */
    mpsort_mpi_newarray(buffer, oldsize, buffer, newsize, packsize, _sort_by_index, sizeof(index[0]), p, comm);

    /* unpack the columns */
    p = buffer;
    for(i = 0; i < newsize; i ++) {
        fastpm_column_set(self, i, (void*) (p + sizeof(index[i])));
        p += packsize;
    }
    fastpm_memory_free(self->mem, buffer);
}
