#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

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

static void
_default_abort(FastPMMemory * m, void * userdata)
{
    fastpm_memory_dump_status(m, 2);
    abort();
}

/*
static void
_default_peak(FastPMMemory * m, void * userdata)
{
    fastpm_memory_dump_status(m, 1);
}
*/
void
fastpm_memory_init(FastPMMemory * m, size_t total_bytes)
{
    int pool;
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        m->pools[pool] = &head;
    }
    m->abortfunc = NULL;
    m->peakfunc = NULL;
    m->userdata = NULL;

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

    m->free_bytes = total_bytes;
    m->total_bytes = total_bytes;
    m->used_bytes = 0;
    m->peak_bytes = 0;
}

void
fastpm_memory_set_handlers(FastPMMemory * m, 
        fastpm_memory_func abortfunc, 
        fastpm_memory_func peakfunc, 
        void * userdata)
{
    m->abortfunc = abortfunc;
    m->peakfunc = peakfunc;

    m->userdata = userdata;
}

static void
_sys_abort(FastPMMemory * m)
{
    if(m->abortfunc) {
        m->abortfunc(m, m->userdata);
    } else {
        _default_abort(m, m->userdata);
    }
}

void
fastpm_memory_tag(FastPMMemory * m, void * p, const char * tag)
{
    MemoryBlock * entry;
    int pool = 0; 
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        for(entry = m->pools[pool]; entry != &head; entry = entry->prev) {
            if(entry->p == p) break;
        }
        if(entry->p == p) {
            /* found in this pool */
            strncpy(entry->tag, tag, 120);
            return;
        }
    }
    /* not found, die */
    if(entry->p != p) _sys_abort(m);
}

void
fastpm_memory_destroy(FastPMMemory * m)
{
    if(m->base0)
        free(m->base0);
    /* check for unrestored pools */
    int pool;
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        if(m->pools[pool] != &head) {
            /* leak !*/
            _sys_abort(m);
        }
    }
}

void
fastpm_memory_dump_status_str(FastPMMemory * m, char * buffer, int n)
{
    /* avoid stdio / malloc. sprintf may still do malloc unfortunately.*/
    MemoryBlock * entry;
    char * buf = buffer;

    n --;
    buf[n] = 0;

    const char P[] = "SHF??????";
    int pool;
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        for(entry = m->pools[pool]; entry != &head; entry = entry->prev) {
            snprintf(buf, n, "%c 0x%016tx : %010td : %s\n", P[pool], (ptrdiff_t) entry->p, entry->size, entry->tag);
            buf += strlen(buf);
            n -= strlen(buf);
            if(n < 0) break;
        }
    }
}
void
fastpm_memory_dump_status(FastPMMemory * m, int fd)
{
    /* avoid stdio / malloc. sprintf may still do malloc unfortunately.*/
    MemoryBlock * entry;
    char buf[1024];

    const char P[] = "SHF??????";
    int pool;
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        for(entry = m->pools[pool]; entry != &head; entry = entry->prev) {
            sprintf(buf, "%c 0x%016tx : %010td : %s\n", P[pool], (ptrdiff_t) entry->p, entry->size, entry->tag);
            write(fd, buf, strlen(buf));
        }
    }
}

static void *
_sys_malloc(FastPMMemory *m, size_t s)
{

    void * p = malloc(s);
    if(p == NULL) {
        _sys_abort(m);
    }
    return p;
}

static void *
fastpm_memory_alloc0(FastPMMemory * m, size_t s, enum FastPMMemoryLocation pool)
{
    s = _align(s, m->alignment);
    if(m->free_bytes <= s) {
        _sys_abort(m);
    }
    m->used_bytes += s;
    m->free_bytes -= s;
    if(m->used_bytes > m->peak_bytes) {
        m->peak_bytes = m->used_bytes;
        if(m->peakfunc)
            m->peakfunc(m, m->userdata);
    }
    MemoryBlock * entry = _sys_malloc(m, sizeof(head));

    int loc = pool;
    if(m->base0 == NULL) { /* allocate from floating but account in the requested loc */
        loc = FASTPM_MEMORY_FLOATING;
    }
    entry->size = s;

    switch(loc) {
        case FASTPM_MEMORY_HEAP:
            entry->p = m->base;
            m->base += s;
        break;
        case FASTPM_MEMORY_STACK:
            entry->p = m->top - s;
            m->top -= s;
        break;
        case FASTPM_MEMORY_FLOATING:
            entry->p = _sys_malloc(m, s);
        break;
    }

    /* add to the pool */
    entry->prev = m->pools[pool];
    m->pools[pool] = entry;

    return entry->p;
}

void *
fastpm_memory_alloc_details(FastPMMemory * m, const char * name,
        size_t s, enum FastPMMemoryLocation pool, const char * file, const int line)
{
    char buf[80];
    sprintf(buf, "%20s: %20s:%d", name, file, line);
    void * r = fastpm_memory_alloc0(m, s, pool);
    fastpm_memory_tag(m, r, buf);
    return r;
}


MemoryBlock * _delist(MemoryBlock ** start, void * p, int * isfirst)
{
    MemoryBlock * entry = NULL;
    MemoryBlock * last= NULL;

    for(entry=*start;
        entry != &head; last=entry, entry=entry->prev) {
        if(entry->p == p) {
            break;
        }
    }
    *isfirst = entry == *start;
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
    int pool;
    for(pool = 0; pool < FASTPM_MEMORY_MAX; pool++) {
        entry = _delist(&m->pools[pool], p, &isfirst);
        if(entry) {
            if(pool != FASTPM_MEMORY_FLOATING && !isfirst) {
                _sys_abort(m);
            }
            int loc = pool;
            /* unbacked */
            if(m->base0 == NULL) {
                loc = FASTPM_MEMORY_FLOATING;
            }
            switch(loc) {
                case FASTPM_MEMORY_STACK:
                    m->top += entry->size;
                break;
                case FASTPM_MEMORY_HEAP:
                    m->base -= entry->size;
                break;
                case FASTPM_MEMORY_FLOATING:
                    free(entry->p);
                break;
            }
            goto exit;
        }
    }
exit:
    m->used_bytes -= entry->size;
    m->free_bytes += entry->size;

    free(entry);

    return;
}
