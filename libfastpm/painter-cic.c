#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"
/* paint and read out */

static double
cic_readout_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], int diffdir);

static void
cic_paint_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight, int diffdir);

void
fastpm_painter_init_cic(FastPMPainter * painter) {
    painter->readout = cic_readout_tuned;
    painter->paint = cic_paint_tuned;
}


static inline double WRtPlus(FastPMFloat * const d, 
        const int i, const int j, const int k, const double f, PM * pm)
{
#pragma omp atomic
    d[k * pm->IRegion.strides[2] + j * pm->IRegion.strides[1] + i * pm->IRegion.strides[0]] += f;
    return f;
}
static inline double REd(FastPMFloat const * const d, const int i, const int j, const int k, const double w, PM * pm)
{
    return d[k * pm->IRegion.strides[2] + j * pm->IRegion.strides[1] + i * pm->IRegion.strides[0]] * w;
}

static void
cic_paint_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight, int diffdir)
{
    PM * pm = painter->pm;
    int d;

    double XYZ[3];
    int IJK[3];
    int IJK1[3];
    double D[3];
    double T[3];

    for(d = 0; d < 3; d ++) {
        XYZ[d] = pos[d]*pm->InvCellSize[d];
        // without floor, -1 < X < 0 is mapped to I=0
        IJK[d] = (int) floor(XYZ[d]);
        IJK1[d] = IJK[d] + 1;
    };

    for(d = 0; d < 3; d ++) {
        D[d] = XYZ[d] - IJK[d];
        T[d] = 1. - D[d];
    }

    if(diffdir >= 0) {
        D[diffdir] = pm->InvCellSize[diffdir];
        T[diffdir] = -pm->InvCellSize[diffdir];
    }


    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    for(d = 0; d < 3; d ++) {
        while(UNLIKELY(IJK[d] < 0)) IJK[d] += pm->Nmesh[d];
        while(UNLIKELY(IJK[d] >= pm->Nmesh[d])) IJK[d] -= pm->Nmesh[d];
        while(UNLIKELY(IJK1[d] < 0)) IJK1[d] += pm->Nmesh[d];
        while(UNLIKELY(IJK1[d] >= pm->Nmesh[d])) IJK1[d] -= pm->Nmesh[d];
    }

    /* start[2] == 0 */
    for(d = 0; d < 2; d ++) {
        IJK[d] -= pm->IRegion.start[d];
        IJK1[d] -= pm->IRegion.start[d];
    }

    D[1] *= weight;
    T[1] *= weight;

    if(LIKELY(0 <= IJK[0] && IJK[0] < pm->IRegion.size[0])) {
        if(LIKELY(0 <= IJK[1] && IJK[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK[0], IJK[1],  IJK[2],  T[2]*T[0]*T[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK[0], IJK[1],  IJK1[2], D[2]*T[0]*T[1], pm);
        }
        if(LIKELY(0 <= IJK1[1] && IJK1[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK[0], IJK1[1], IJK[2],  T[2]*T[0]*D[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK[0], IJK1[1], IJK1[2], D[2]*T[0]*D[1], pm);
        }
    }

    if(LIKELY(0 <= IJK1[0] && IJK1[0] < pm->IRegion.size[0])) {
        if(LIKELY(0 <= IJK[1] && IJK[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK1[0], IJK[1],  IJK[2],  T[2]*D[0]*T[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK1[0], IJK[1],  IJK1[2], D[2]*D[0]*T[1], pm);
        }
        if(LIKELY(0 <= IJK1[1] && IJK1[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK1[0], IJK1[1], IJK[2],  T[2]*D[0]*D[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                WRtPlus(canvas, IJK1[0], IJK1[1], IJK1[2], D[2]*D[0]*D[1], pm);
        }
    }
}

static double
cic_readout_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], int diffdir)
{
    PM * pm = painter->pm;

    int d;

    double XYZ[3];
    int IJK[3];
    int IJK1[3];
    double D[3];
    double T[3];

    for(d = 0; d < 3; d ++) {
        XYZ[d] = pos[d]*pm->InvCellSize[d];
        // without floor, -1 < X < 0 is mapped to I=0
        IJK[d] = (int) floor(XYZ[d]);
        IJK1[d] = IJK[d] + 1;
    };

    for(d = 0; d < 3; d ++) {
        D[d] = XYZ[d] - IJK[d];
        T[d] = 1. - D[d];
    }

    if(diffdir >= 0) {
        D[diffdir] = pm->InvCellSize[diffdir];
        T[diffdir] = -pm->InvCellSize[diffdir];
    }


    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    for(d = 0; d < 3; d ++) {
        while(UNLIKELY(IJK[d] < 0)) IJK[d] += pm->Nmesh[d];
        while(UNLIKELY(IJK[d] >= pm->Nmesh[d])) IJK[d] -= pm->Nmesh[d];
        while(UNLIKELY(IJK1[d] < 0)) IJK1[d] += pm->Nmesh[d];
        while(UNLIKELY(IJK1[d] >= pm->Nmesh[d])) IJK1[d] -= pm->Nmesh[d];
    }

    /* start[2] == 0 */
    for(d = 0; d < 2; d ++) {
        IJK[d] -= pm->IRegion.start[d];
        IJK1[d] -= pm->IRegion.start[d];
    }

    double value = 0;

    if(LIKELY(0 <= IJK[0] && IJK[0] < pm->IRegion.size[0])) {
        if(LIKELY(0 <= IJK[1] && IJK[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK[0], IJK[1],  IJK[2],  T[2]*T[0]*T[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK[0], IJK[1],  IJK1[2], D[2]*T[0]*T[1], pm);
        }
        if(LIKELY(0 <= IJK1[1] && IJK1[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK[0], IJK1[1], IJK[2],  T[2]*T[0]*D[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK[0], IJK1[1], IJK1[2], D[2]*T[0]*D[1], pm);
        }
    }

    if(LIKELY(0 <= IJK1[0] && IJK1[0] < pm->IRegion.size[0])) {
        if(LIKELY(0 <= IJK[1] && IJK[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK1[0], IJK[1],  IJK[2],  T[2]*D[0]*T[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK1[0], IJK[1],  IJK1[2], D[2]*D[0]*T[1], pm);
        }
        if(LIKELY(0 <= IJK1[1] && IJK1[1] < pm->IRegion.size[1])) {
            if(LIKELY(0 <= IJK[2] && IJK[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK1[0], IJK1[1], IJK[2],  T[2]*D[0]*D[1], pm);
            if(LIKELY(0 <= IJK1[2] && IJK1[2] < pm->IRegion.size[2]))
                value += REd(canvas, IJK1[0], IJK1[1], IJK1[2], D[2]*D[0]*D[1], pm);
        }
    }
    return value;
}

