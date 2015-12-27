#include <stddef.h>
#include <stdint.h>

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

FASTPM_END_DECLS

#include "solver-2lpt.h"
#include "solver-pm.h"
#include "utils.h"
#include "utils-ps.h"
#include "utils-ic.h"

#include "pmapi.h"
