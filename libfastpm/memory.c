#include <stdlib.h>
#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/memory.h>
#include <fastpm/logging.h>


struct MemoryBlock {
    void * p;
    size_t size;
    MemoryBlock * prev; /* pointer to previous block */
    char tag[128]; /* tag */
};

static MemoryBlock head = { .p = NULL, .size = 0, .prev = NULL};

static size_t
_align(size_t old, size_t alignment)
{
    size_t n;
    if(alignment == 0) return old;

    n = (old / alignment) * alignment;

    while (n <= old) n += alignment;
    return n;
}

void
fastpm_memory_init(FastPMMemory * m, size_t total_bytes, int allow_unordered)
{
    m->heap = &head;
    m->stack = &head;

    if(total_bytes > 0) {
        total_bytes = _align(total_bytes, m->alignment);

        m->base0 = (char*) malloc(total_bytes + m->alignment);
        m->base = (char*) _align((size_t)(m->base0), m->alignment);
        m->top = m->base + total_bytes;
    } else {
        total_bytes = 0x8000000000000000;
        m->base0 = NULL;
        m->base = NULL;
        m->top = NULL;

    }

    m->allow_unordered = 1;
    m->free_bytes = total_bytes;
    m->total_bytes = total_bytes;
    m->used_bytes = 0;
    m->peak_bytes = 0;
}

void
fastpm_memory_tag(FastPMMemory * m, void * p, const char * tag)
{
    MemoryBlock * entry;
    for(entry = m->stack; entry != &head; entry = entry->prev) {
        if(entry->p == p) break;
    }
    for(entry = m->heap; entry != &head; entry = entry->prev) {
        if(entry->p == p) break;
    }
    if(entry->p != p) abort();

    strncpy(entry->tag, tag, 120);
}

void
fastpm_memory_destroy(FastPMMemory * m)
{
    if(m->base0)
        free(m->base0);
    if(m->stack != &head) {
        /* leak !*/
        abort();
    }
    if(m->heap != &head) {
        /* leak !*/
        abort();
    }
}

void *
fastpm_memory_alloc(FastPMMemory * m, size_t s, enum FastPMMemoryLocation loc)
{
    s = _align(s, m->alignment);
    if(m->free_bytes <= s) {
        abort();
    }
    m->used_bytes += s;
    m->free_bytes -= s;
    if(m->used_bytes > m->peak_bytes) {
        m->peak_bytes = m->used_bytes;
    }
    MemoryBlock * entry = malloc(sizeof(head));
    switch(loc) {
        case FASTPM_MEMORY_HEAP:
            if(m->base) {
                entry->p = m->base;
                m->base += s;
            } else {
                entry->p = malloc(s);
            }
            entry->size = s;
            entry->prev = m->heap;
            m->heap = entry;
        break;
        case FASTPM_MEMORY_STACK:
            if(m->top) {
                entry->p = m->top - s;
                m->top -= s;
            } else {
                entry->p = malloc(s);
            }
            entry->size = s;
            entry->prev = m->stack;
            m->stack = entry;
        break;
    }
    return entry->p;
}

MemoryBlock * _delist(MemoryBlock ** start, void * p, int * isfirst)
{
    MemoryBlock * entry = NULL;
    MemoryBlock * last= NULL;
    *isfirst = entry == *start;

    for(entry=*start;
        entry != &head; last=entry, entry=entry->prev) {
        if(entry->p == p) {
            break;
        }
    }
    if(entry == &head) return NULL;

    if(last) {
        last->prev = entry->prev;
    } else {
        *start = entry->prev;
    }
    return entry;
}


void
fastpm_memory_free(FastPMMemory * m, void * p)
{
    MemoryBlock * entry;
    int isfirst = 0;
    entry = _delist(&m->stack, p, &isfirst);
    if(entry) {
        if(!m->allow_unordered && !isfirst) {
            abort();
        }
        if(m->top) {
            m->top += entry->size;
        } else {
            free(entry->p);
        }
        goto exit;
    }

    entry = _delist(&m->heap, p, &isfirst);
    if(entry) {
        if(!m->allow_unordered && !isfirst) {
            abort();
        }
        if(m->base) {
            m->base -= entry->size;
        } else {
            free(entry->p);
        }
        goto exit;
    }

exit:
    m->used_bytes -= entry->size;
    m->free_bytes += entry->size;

    free(entry);

    return;
}
