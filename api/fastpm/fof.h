#ifndef _FASTPM_FOF_H
#define _FASTPM_FOF_H

#define FASTPM_EVENT_HALO "HALO"
typedef struct FastPMFOFFinderPrivate FastPMFOFFinderPrivate;

typedef struct {
    int nmin;
    int periodic;
    int kdtree_thresh;

    FastPMStore * p;
    PM * pm;

    FastPMFOFFinderPrivate * priv;

} FastPMFOFFinder;

/* linking length is in absolute units (aka not b, but b * mean_sep)*/
void
fastpm_fof_init(FastPMFOFFinder * finder,
                double max_linkinglength,
                FastPMStore * store, PM * pm);

/* create a halo catalog from the heap.
 * halos->name shall be set before this.
 *
 * if active is not NULL, only run fof on particles active[i] != 0.
 * returns ihalo to the halo id of each particle in p,
 * halos contains every halo that spans to the rank. primary halos
 * has mask[ihalo] == 1. use subsample to remove non-primary halos.
 *
 * halos is created on the heap, ownership transferred to the caller.
 * ihalo is created ho the stack, ownership transferred to the caller.
 * */
ptrdiff_t *
fastpm_fof_execute(FastPMFOFFinder * finder,
                   double linkinglength,
                   FastPMStore * halos,
                   FastPMParticleMaskType * active);

void
fastpm_fof_subsample_and_relabel(FastPMFOFFinder * finder,
    FastPMStore * halos,
    FastPMParticleMaskType * mask,
    ptrdiff_t * head);

void
fastpm_fof_destroy(FastPMFOFFinder * finder);

void
fastpm_fof_allocate_halos(FastPMStore * halos,
    size_t nhalos,
    FastPMStore * p,
    int include_q,
    MPI_Comm comm);

#endif
