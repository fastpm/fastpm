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

struct FastPMMemory {
    size_t alignment;
    size_t total_bytes;

    MemoryBlock * heap;
    MemoryBlock * stack;
    size_t peak_bytes;
    size_t used_bytes;
    size_t free_bytes;
    char * base0;
    char * top;
    char * base;
};

enum FastPMMemoryLocation {
    FASTPM_MEMORY_HEAP,
    FASTPM_MEMORY_STACK,
};

void
fastpm_memory_init(FastPMMemory * m, size_t total_bytes);

void
fastpm_memory_tag(FastPMMemory * m, void * p, const char * tag);

void
fastpm_memory_destroy(FastPMMemory * m);

void
fastpm_memory_free(FastPMMemory * m, void * p);

void *
fastpm_memory_alloc(FastPMMemory * m, size_t s, enum FastPMMemoryLocation loc);

FASTPM_END_DECLS

#endif
