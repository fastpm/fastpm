#ifndef _FASTPM_FOF_H
#define _FASTPM_FOF_H
typedef struct FastPMFOFFinderPrivate FastPMFOFFinderPrivate;

typedef struct {
    int nmin;
    double linkinglength; /* absolute */
    int periodic;
    FastPMStore halos[1];
    int kdtree_thresh;

    FastPMStore * p;
    PM * pm;

    /* private */
    uint64_t * label;

    FastPMFOFFinderPrivate * priv;

} FastPMFOFFinder;

void
fastpm_fof_init(FastPMFOFFinder * finder, FastPMStore * store, PM * pm);

void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos);

void
fastpm_fof_destroy(FastPMFOFFinder * finder);
#endif
