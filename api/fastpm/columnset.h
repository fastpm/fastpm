#ifndef _FASTPM_COLUMNSET_H_
#define _FASTPM_COLUMNSET_H_

#include <stddef.h>
#include <stdint.h>

/* Subclass of fastpm_column. */
FASTPM_BEGIN_DECLS
typedef struct {
    FastPMColumn base;
    FastPMColumn * columns[20];
    int ncolumns;
} FastPMColumnSet;

void
fastpm_column_set_init(FastPMColumnSet * self,
        FastPMColumn ** columns, int ncolumns
        );

FASTPM_END_DECLS

#endif
