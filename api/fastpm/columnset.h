#ifndef _FASTPM_COLUMNSET_H_
#define _FASTPM_COLUMNSET_H_

#include <stddef.h>
#include <stdint.h>

/* Subclass of fastpm_column. */
FASTPM_BEGIN_DECLS
typedef struct {
    FastPMColumn base;
    FastPMColumn * columns[80];
} FastPMColumnSet;

/* create a columnset from NULL terminated list of columns;
 * at most 80 is supported*/

void
fastpm_columnset_init(FastPMColumnSet * self, FastPMColumn ** columns);

FASTPM_END_DECLS

#endif
