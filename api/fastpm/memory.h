#ifndef __FASTPM_MEMORY_H__
#define __FASTPM_MEMORY_H__

#ifndef FASTPM_BEGIN_DECLS
#ifdef __cplusplus
#define FASTPM_BEGIN_DECLS extern "C" {
#define FASTPM_END_DECLS }
#else
#define FASTPM_BEGIN_DECLS
#define FASTPM_END_DECLS
#endif
#endif

FASTPM_BEGIN_DECLS

typedef struct MemoryBlock MemoryBlock;
typedef struct FastPMMemory FastPMMemory;
typedef void (*fastpm_memory_func)(FastPMMemory * m, void * userdata);

struct FastPMMemory {
    size_t alignment;
    size_t total_bytes;

    MemoryBlock * pools[8];
    size_t peak_bytes;
    size_t used_bytes;
    size_t free_bytes;
    char * base0;
    char * top;
    char * base;

    void * userdata;
    fastpm_memory_func abortfunc;
    fastpm_memory_func peakfunc;
};

enum FastPMMemoryLocation {
    FASTPM_MEMORY_HEAP,
    FASTPM_MEMORY_STACK,
    FASTPM_MEMORY_FLOATING,
    FASTPM_MEMORY_MAX,
};

void
fastpm_memory_init(FastPMMemory * m, size_t total_bytes);

void
fastpm_memory_set_handlers(FastPMMemory * m, fastpm_memory_func abortfunc, fastpm_memory_func peakfunc, void * userdata);

void
fastpm_memory_tag(FastPMMemory * m, void * p, const char * tag);

void
fastpm_memory_destroy(FastPMMemory * m);

void
fastpm_memory_free(FastPMMemory * m, void * p);

void *
fastpm_memory_alloc_details(FastPMMemory * m, const char * name, size_t s, enum FastPMMemoryLocation loc, const char * file, const int line);

void
fastpm_memory_dump_status(FastPMMemory * m, int fd);

void
fastpm_memory_dump_status_str(FastPMMemory * m, char * buf, int n);

#define fastpm_memory_alloc(m, name, s, loc) fastpm_memory_alloc_details(m, name, s, loc, __FILE__, __LINE__)

FASTPM_END_DECLS

#endif
