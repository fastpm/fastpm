#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include "pmpfft.h"

/* paint and read out */

static inline double WRtPlus(real_t * const d, 
        const int i, const int j, const int k, const double f, PM * pm)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
    d[k * pm->IRegion.strides[2] + j * pm->IRegion.strides[1] + i * pm->IRegion.strides[0]] += f;
    return f;
}
static inline double REd(real_t const * const d, const int i, const int j, const int k, const double w, PM * pm)
{
    return d[k * pm->IRegion.strides[2] + j * pm->IRegion.strides[1] + i * pm->IRegion.strides[0]] * w;
}

static inline void pm_paint_pos_tuned(PM * pm, double pos[3], double weight) {
    double X=pos[0]*pm->InvCellSize[0];
    double Y=pos[1]*pm->InvCellSize[1];
    double Z=pos[2]*pm->InvCellSize[2];

    int I=(int) floor(X); // without floor, -1 < X < 0 is mapped to I=0
    int J=(int) floor(Y);          // Assumes Y,Z are positive
    int K=(int) floor(Z);
    double D1=X-((double) I);
    double D2=Y-((double) J);
    double D3=Z-((double) K);
    double T1=1.-D1;
    double T2=1.-D2;
    double T3=1.-D3;

    double T2W =T2*weight;
    double D2W =D2*weight;

    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    while(UNLIKELY(I >= pm->Nmesh[0])) I -= pm->Nmesh[0];
    while(UNLIKELY(J >= pm->Nmesh[1])) J -= pm->Nmesh[1];
    while(UNLIKELY(K >= pm->Nmesh[2])) K -= pm->Nmesh[2];

    int I1=I+1; while(UNLIKELY(I1 >= pm->Nmesh[0])) I1-=pm->Nmesh[0];
    int J1=J+1; while(UNLIKELY(J1 >= pm->Nmesh[1])) J1-=pm->Nmesh[1]; // assumes y,z < BoxSize
    int K1=K+1; while(UNLIKELY(K1 >= pm->Nmesh[2])) K1-=pm->Nmesh[2];

    I -= pm->IRegion.start[0];
    I1 -= pm->IRegion.start[0];
    J -= pm->IRegion.start[0];
    J1 -= pm->IRegion.start[0];

    if(LIKELY(0 <= I && I < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            WRtPlus(pm->workspace, I, J,  K,  T3*T1*T2W, pm);
            WRtPlus(pm->workspace, I, J,  K1, D3*T1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            WRtPlus(pm->workspace, I, J1, K,  T3*T1*D2W, pm);
            WRtPlus(pm->workspace, I, J1, K1, D3*T1*D2W, pm);
        }
    }

    if(LIKELY(0 <= I1 && I1 < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            WRtPlus(pm->workspace, I1, J,  K,  T3*D1*T2W, pm);
            WRtPlus(pm->workspace, I1, J,  K1, D3*D1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            WRtPlus(pm->workspace, I1, J1, K,  T3*D1*D2W, pm);
            WRtPlus(pm->workspace, I1, J1, K1, D3*D1*D2W, pm);
        }
    }
}

static inline double pm_readout_pos_tuned(PM * pm, double pos[3]) {
    double X=pos[0]*pm->InvCellSize[0];
    double Y=pos[1]*pm->InvCellSize[1];
    double Z=pos[2]*pm->InvCellSize[2];

    int I=(int) floor(X); // without floor, -1 < X < 0 is mapped to I=0
    int J=(int) floor(Y);          // Assumes Y,Z are positive
    int K=(int) floor(Z);
    double D1=X-((double) I);
    double D2=Y-((double) J);
    double D3=Z-((double) K);
    double T1=1.-D1;
    double T2=1.-D2;
    double T3=1.-D3;

    double T2W =T2;
    double D2W =D2;

    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    while(UNLIKELY(I >= pm->Nmesh[0])) I -= pm->Nmesh[0];
    while(UNLIKELY(J >= pm->Nmesh[1])) J -= pm->Nmesh[1];
    while(UNLIKELY(K >= pm->Nmesh[2])) K -= pm->Nmesh[2];

    int I1=I+1; while(UNLIKELY(I1 >= pm->Nmesh[0])) I1-=pm->Nmesh[0];
    int J1=J+1; while(UNLIKELY(J1 >= pm->Nmesh[1])) J1-=pm->Nmesh[1]; // assumes y,z < BoxSize
    int K1=K+1; while(UNLIKELY(K1 >= pm->Nmesh[2])) K1-=pm->Nmesh[2];

    I -= pm->IRegion.start[0];
    I1 -= pm->IRegion.start[0];
    J -= pm->IRegion.start[0];
    J1 -= pm->IRegion.start[0];

    double value = 0;

    if(LIKELY(0 <= I && I < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            value += REd(pm->workspace, I, J,  K,  T3*T1*T2W, pm);
            value += REd(pm->workspace, I, J,  K1, D3*T1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            value += REd(pm->workspace, I, J1, K,  T3*T1*D2W, pm);
            value += REd(pm->workspace, I, J1, K1, D3*T1*D2W, pm);
        }
    }

    if(LIKELY(0 <= I1 && I1 < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            value += REd(pm->workspace, I1, J,  K,  T3*D1*T2W, pm);
            value += REd(pm->workspace, I1, J,  K1, D3*D1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            value += REd(pm->workspace, I1, J1, K,  T3*D1*D2W, pm);
            value += REd(pm->workspace, I1, J1, K1, D3*D1*D2W, pm);
        }
    }
    return value;
}

static inline void pm_paint_pos(PM * pm, double pos[3], double weight) {
    /* to pm->workspace */
    int n;
    double gpos[3];
    int ipos[3];
    float k[3][2];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->InvCellSize[d];
        ipos[d] = floor(gpos[d]);
        k[d][0] = 1 + ipos[d] - gpos[d];
        k[d][1] = gpos[d] - ipos[d];
        ipos[d] -= pm->IRegion.start[d];
    }
    for(n = 0; n < 8; n ++) {
        float kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n>>d) & 1;
            int targetpos = ipos[d] + rel;
            kernel *= k[d][rel];
            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            ind += pm->IRegion.strides[d] * targetpos;
            if(targetpos < 0 || targetpos >= pm->IRegion.size[d]) 
                goto outside;
        }
        pm->workspace[ind] += weight * kernel;

    outside:
        continue;
    }
    return;
}

static double pm_readout_pos(PM * pm, double pos[3]) {
    /* from pm->workspace */
    double value = 0;
    int n;
    double gpos[3];
    int ipos[3];
    float k[3][2];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->InvCellSize[d];
        ipos[d] = floor(gpos[d]);
        k[d][0] = 1 + ipos[d] - gpos[d];
        k[d][1] = gpos[d] - ipos[d];
    }
    for(n = 0; n < 8; n ++) {
        float kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n>>d) & 1;

            kernel *= k[d][rel];

            int targetpos = ipos[d] + rel - pm->IRegion.start[d];

            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) {
                goto outside;
            } 
            ind += pm->IRegion.strides[d] * targetpos;
        }
        value += kernel * pm->workspace[ind];
outside:
        continue;
    }
    return value;
}

void pm_paint(PM * pm, void * pdata, ptrdiff_t size) {
    ptrdiff_t i;
    memset(pm->workspace, 0, sizeof(pm->workspace[0]) * pm->allocsize);
    for (i = 0; i < size; i ++) {
        double pos[3];
        pm->iface.get_position(pdata, i, pos);
        pm_paint_pos_tuned(pm, pos, 1.0);
//        pm_paint_pos(pm, pos, 1.0);
    }
}

double pm_readout_one(PM * pm, void * pdata, ptrdiff_t i) {
    double pos[3];
    pm->iface.get_position(pdata, i, pos);
 //   return pm_readout_pos(pm, pos);
    return pm_readout_pos_tuned(pm, pos);
}
