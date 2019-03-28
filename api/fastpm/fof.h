#ifndef _FASTPM_FOF_H
#define _FASTPM_FOF_H

#define FASTPM_EVENT_HALO "HALO"
typedef struct FastPMFOFFinderPrivate FastPMFOFFinderPrivate;

typedef struct {
    FastPMEvent base;
    FastPMStore * halos;
    FastPMStore * p;
    ptrdiff_t * ihalo; /* halo index in halo per particle. */
} FastPMHaloEvent;

typedef struct {
    int nmin;
    double linkinglength; /* absolute */
    int periodic;
    int kdtree_thresh;

    FastPMStore * p;
    PM * pm;

    /* private */
    uint64_t * label;

    FastPMFOFFinderPrivate * priv;

    FastPMEventHandler * event_handlers;
} FastPMFOFFinder;

void
fastpm_fof_init(FastPMFOFFinder * finder, FastPMStore * store, PM * pm);

/* create a halo catalog from the heap. halos->name shall be set before this */
void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos);

void
fastpm_fof_destroy(FastPMFOFFinder * finder);
#endif
