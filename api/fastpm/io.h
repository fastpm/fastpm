#ifndef __FASTPM_IO_H__
#define __FASTPM_IO_H__
#include <bigfile.h>
#include <fastpm/histogram.h>

FASTPM_BEGIN_DECLS

typedef void (*FastPMSnapshotSorter)(const void * ptr, void * radix, void * arg);

void
FastPMSnapshotSortByID(const void * ptr, void * radix, void * arg);

void
FastPMSnapshotSortByLength(const void * ptr, void * radix, void * arg);

void
FastPMSnapshotSortByAEmit(const void * ptr, void * radix, void * arg);

void
fastpm_sort_snapshot(FastPMStore * p, MPI_Comm comm, FastPMSnapshotSorter sorter, int redistribute);

int
fastpm_store_write(FastPMStore * p,
        const char * filebase,
        const char * mode,
        int Nwriters,
        MPI_Comm comm
);

int
fastpm_store_read(FastPMStore * p,
        const char * filebase,
        int Nwriters,
        MPI_Comm comm
);

void
write_snapshot_header(FastPMSolver * fastpm,
    const char * filebase, MPI_Comm comm);

void
read_snapshot_header(FastPMSolver * fastpm,
    const char * filebase, double * aout, MPI_Comm comm);

void
write_snapshot_attr(const char * filebase,
    const char * dataset,
    const char * attrname,
    void * buf,
    const char * dtype,
    size_t nmemb,
    MPI_Comm comm);

int
read_snapshot(FastPMSolver * fastpm, FastPMStore * p, const char * filebase);

int
write_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters);

int
read_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters);

size_t
read_angular_grid(FastPMStore * store,
        const char * filename,
        const double * r,
        const double * aemit,
        const size_t Nr,
        int sampling_factor,
        MPI_Comm comm);

void
write_aemit_hist(const char * fn, const char * ds,
            FastPMHistogram * hist,
            MPI_Comm comm);

int
read_funck(FastPMFuncK * fk, const char filename[], MPI_Comm comm);

/* paints i-th partile in store to i-th particle in map.*/
typedef void (*fastpm_store_hp_paintfunc)(
        FastPMStore * store,
        ptrdiff_t i,
        FastPMStore * map,
        void * userdata);

void
fastpm_snapshot_paint_hpmap(FastPMStore * p,
        int64_t nside,
        int64_t nslice,
        fastpm_store_hp_paintfunc paintfunc,
        void * userdata,
        FastPMStore * map,
        MPI_Comm comm
);

void
fastpm_store_histogram_aemit_sorted(FastPMStore * store,
        FastPMHistogram * hist,
        MPI_Comm comm);

FASTPM_END_DECLS
#endif
