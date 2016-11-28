#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

static void
_fastpm_column_set_destroy(FastPMColumn * base);

static void
_fastpm_column_set_set(FastPMColumn * base, ptrdiff_t i, void * src);

static void
_fastpm_column_set_get(FastPMColumn * base, ptrdiff_t i, void * dest);

void
fastpm_column_set_init(FastPMColumnSet * self,
        FastPMColumn ** columns, int ncolumns
        )
{
    FastPMColumn * base = &self->base;
    base->mem = columns[0]->mem;
    int i;
    base->nmemb = 1;
    base->size = columns[0]->size;
    base->maxsize = columns[0]->maxsize;
    base->elsize = 0;
    base->buffer = NULL;
    self->ncolumns = ncolumns;

    for(i = 0; i < ncolumns; i ++) {
        base->elsize += columns[i]->elsize * columns[i]->nmemb;
        if(columns[i]->maxsize != base->maxsize) {
            fastpm_raise(-1, "mismatched maxsize (%td != %td)\n", columns[i]->maxsize, base->maxsize);
        }
        if(columns[i]->size != base->size) {
            fastpm_raise(-1, "mismatched size (%td != %td)\n", columns[i]->size, base->size);
        }
        self->columns[i] = columns[i];
    }
    base->to_double = NULL;
    base->from_double = NULL;
    base->get = _fastpm_column_set_get;
    base->set = _fastpm_column_set_set;
    base->destroy = _fastpm_column_set_destroy;
}

static void
_fastpm_column_set_get(FastPMColumn * base, ptrdiff_t i, void * dest)
{
    int c;
    FastPMColumnSet * self = (FastPMColumnSet*) base;
    char * buf = (char*) dest;
    for(c = 0; c < self->ncolumns; c++) {
        fastpm_column_get(self->columns[c], i, buf);
        buf += self->columns[c]->elsize * self->columns[c]->nmemb;
    }
}

static void
_fastpm_column_set_set(FastPMColumn * base, ptrdiff_t i, void * src)
{
    int c; 
    FastPMColumnSet * self = (FastPMColumnSet*) base;
    char * buf = (char*) src;
    for(c = 0; c < self->ncolumns; c++) {
        fastpm_column_set(self->columns[c], i, buf);
        buf += self->columns[c]->elsize * self->columns[c]->nmemb;
    }
}

static void
_fastpm_column_set_destroy(FastPMColumn * base)
{
    return;
}
