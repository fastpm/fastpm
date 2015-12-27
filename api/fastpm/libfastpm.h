#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PM PM;
typedef struct PMStore PMStore;

typedef struct VPM VPM;
typedef struct PMGhostData PMGhostData;

typedef struct {
    void * (*malloc )(size_t);
    void   (*free   )(void *);
    void   (*get_position)(void * pdata, ptrdiff_t index, double pos[3]);
    size_t (*pack)  (void * pdata, ptrdiff_t index, void * packed, int attributes);
    void   (*unpack)(void * pdata, ptrdiff_t index, void * packed, int attributes);
    void   (*reduce)(void * pdata, ptrdiff_t index, void * packed, int method);
} PMIFace;

typedef struct {
    /* in units of real numbers, not bytes. */
    ptrdiff_t start[3];
    ptrdiff_t size[3];
    ptrdiff_t strides[3]; 
    ptrdiff_t total;
} PMRegion;

struct PMStore {
    PMIFace iface;

    int attributes; /* bit flags of allocated attributes */

    double (* x)[3];
    float (* q)[3];
    float (* v)[3];
    float (* acc)[3];
    float (* dx1)[3];
    float (* dx2)[3];
    uint64_t * id;
    size_t np;
    size_t np_upper;
    double a_x;
    double a_v;
};

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

#include "2lpt.h"
#include "pm.h"
#include "utils.h"

#include "pmkiter.h"
#include "pmpowerspectrum.h"

#ifdef __cplusplus
}
#endif
