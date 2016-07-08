#ifndef __FASTPM_H__
#define __FASTPM_H__

#include <stddef.h>
#include <stdint.h>
#include <mpi.h>

#ifndef FASTPM_BEGIN_DECLS
#ifdef __cplusplus
#define FASTPM_BEGIN_DECLS extern "C" {
#define FASTPM_END_DECLS }
#else
#define FASTPM_BEGIN_DECLS
#define FASTPM_END_DECLS
#endif
#endif

FASTPM_BEGIN_DECLS

typedef struct PM PM;
typedef struct PMStore PMStore;
typedef struct FastPMPainter FastPMPainter;

#ifndef FASTPM_FFT_PRECISION
#define FASTPM_FFT_PRECISION 32
#endif

#if FASTPM_FFT_PRECISION == 64
    typedef double FastPMFloat;
#elif FASTPM_FFT_PRECISION == 32
    typedef float FastPMFloat;
#else
    #error FASTPM_FFT_PRECISION must be 32 or 64
#endif

void libfastpm_init();
void libfastpm_cleanup();
void libfastpm_set_memory_bound(size_t size);

FASTPM_END_DECLS

typedef double (*fastpm_fkfunc)(double k, void * data);
typedef void   (*fastpm_posfunc)(PMStore * p, ptrdiff_t index, double pos[3]);
typedef double (*fastpm_kernelfunc)(double x, int support);
#define fastpm_pkfunc fastpm_fkfunc

#include "painter.h"
#include "memory.h"
#include "pmapi.h"

#include "transfer.h"
#include "utils.h"
#include "initialcondition.h"
#include "pngaussian.h"
#include "powerspectrum.h"
#include "solver-pm.h"

/* Following functions are internal API */

FASTPM_BEGIN_DECLS

FastPMMemory * _libfastpm_get_gmem();

FASTPM_END_DECLS

#endif
