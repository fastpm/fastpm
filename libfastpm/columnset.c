#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

static void
_fastpm_columnset_resize(FastPMColumn * base, size_t newsize);

static void
_fastpm_columnset_destroy(FastPMColumn * base);

static void
_fastpm_columnset_set(FastPMColumn * base, ptrdiff_t i, void * src);

static void
_fastpm_columnset_get(FastPMColumn * base, ptrdiff_t i, void * dest);

void
fastpm_columnset_init(FastPMColumnSet * self, FastPMColumn ** columns)
{
    FastPMColumn * base = &self->base;
    base->mem = columns[0]->mem;
    base->nmemb = 1;
    base->size = columns[0]->size;
    base->maxsize = columns[0]->maxsize;
    base->elsize = 0;
    base->buffer = NULL;

    int c;
    FastPMColumn * column;
    for(c = 0; NULL != (column = columns[c]); c++) {
        base->elsize += column->rowsize;
        if(column->maxsize < base->maxsize) {
            base->maxsize = column->maxsize;
        }
        if(column->size != base->size) {
            fastpm_raise(-1, "mismatched size (%td != %td)\n", column->size, base->size);
        }
        if(c >= 80) {
            fastpm_raise(-1, "fix the code to support more than 80 columns.\n");
        }
        self->columns[c] = column;
    }
    if(base->maxsize < base->size) {
        fastpm_raise(-1, "mismatched maxsize (%td < %td)\n", base->maxsize, base->size);
    }
    base->rowsize = base->elsize;

    /* terminate the list with NULL */
    self->columns[c] = NULL;

    base->to_double = NULL;
    base->from_double = NULL;
    base->get = _fastpm_columnset_get;
    base->set = _fastpm_columnset_set;
    base->destroy = _fastpm_columnset_destroy;
    base->resize = _fastpm_columnset_resize;
}

static void
_fastpm_columnset_resize(FastPMColumn * base, size_t newsize)
{
    base->size = newsize;
    FastPMColumnSet * self = (FastPMColumnSet *) base;
    int c;
    FastPMColumn * column;
    for(c = 0; NULL != (column = self->columns[c]); c++) {
        fastpm_column_resize(column, newsize);
    }

}

static void
_fastpm_columnset_get(FastPMColumn * base, ptrdiff_t i, void * dest)
{
    FastPMColumnSet * self = (FastPMColumnSet*) base;
    char * buf = (char*) dest;
    int c;
    FastPMColumn * column;
    for(c = 0; NULL != (column = self->columns[c]); c++) {
        /* if the column doesn't use any storage (const), skip it */
        if(column->rowsize > 0)
            fastpm_column_get(column, i, buf);
        buf += column->rowsize;
    }
}

static void
_fastpm_columnset_set(FastPMColumn * base, ptrdiff_t i, void * src)
{
    FastPMColumnSet * self = (FastPMColumnSet*) base;
    char * buf = (char*) src;
    int c;
    FastPMColumn * column;
    for(c = 0; NULL != (column = self->columns[c]); c++) {
        /* if the column doesn't use any storage (const), skip it */
        if(column->rowsize > 0)
            fastpm_column_set(column, i, buf);
        buf += column->rowsize;
    }
}

static void
_fastpm_columnset_destroy(FastPMColumn * base)
{
    return;
}


