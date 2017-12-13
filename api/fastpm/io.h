#ifndef __FASTPM_IO_H__
#define __FASTPM_IO_H__
FASTPM_BEGIN_DECLS

typedef void (*FastPMSnapshotSorter)(const void * ptr, void * radix, void * arg);

void
FastPMSnapshotSortByID(const void * ptr, void * radix, void * arg);

void
FastPMSnapshotSortByAEmit(const void * ptr, void * radix, void * arg);

int 
write_snapshot(FastPMSolver * fastpm,
        FastPMStore * p,
        char * filebase,
        char * parameters,
        int Nwriters,
        FastPMSnapshotSorter sorter;
    );

int 
read_snapshot(FastPMSolver * fastpm, FastPMStore * p, char * filebase);

int
write_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters);

int
read_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters);

FASTPM_END_DECLS
#endif
