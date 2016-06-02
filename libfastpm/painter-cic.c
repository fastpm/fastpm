#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"
/* paint and read out */

static double
cic_readout_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3]);

static void
cic_paint_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight);

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
cic_paint_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight)
{
    PM * pm = painter->pm;
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
    double T3=1.-D3;

    double D2W = D2*weight;
    double T2W = weight -D2W;

    int I1=I+1; 
    int J1=J+1; 
    int K1=K+1; 

    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    while(UNLIKELY(I < 0)) I += pm->Nmesh[0];
    while(UNLIKELY(J < 0)) J += pm->Nmesh[1];
    while(UNLIKELY(K < 0)) K += pm->Nmesh[2];
    while(UNLIKELY(I >= pm->Nmesh[0])) I -= pm->Nmesh[0];
    while(UNLIKELY(J >= pm->Nmesh[1])) J -= pm->Nmesh[1];
    while(UNLIKELY(K >= pm->Nmesh[2])) K -= pm->Nmesh[2];


    while(UNLIKELY(I1 < 0)) I1 += pm->Nmesh[0];
    while(UNLIKELY(J1 < 0)) J1 += pm->Nmesh[1];
    while(UNLIKELY(K1 < 0)) K1 += pm->Nmesh[2];
    while(UNLIKELY(I1 >= pm->Nmesh[0])) I1-=pm->Nmesh[0];
    while(UNLIKELY(J1 >= pm->Nmesh[1])) J1-=pm->Nmesh[1];
    while(UNLIKELY(K1 >= pm->Nmesh[2])) K1-=pm->Nmesh[2];

    I -= pm->IRegion.start[0];
    J -= pm->IRegion.start[1];
    I1 -= pm->IRegion.start[0];
    J1 -= pm->IRegion.start[1];

    if(LIKELY(0 <= I && I < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                WRtPlus(canvas, I, J,  K,  T3*T1*T2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                WRtPlus(canvas, I, J,  K1, D3*T1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                WRtPlus(canvas, I, J1, K,  T3*T1*D2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                WRtPlus(canvas, I, J1, K1, D3*T1*D2W, pm);
        }
    }

    if(LIKELY(0 <= I1 && I1 < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                WRtPlus(canvas, I1, J,  K,  T3*D1*T2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                WRtPlus(canvas, I1, J,  K1, D3*D1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                WRtPlus(canvas, I1, J1, K,  T3*D1*D2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                WRtPlus(canvas, I1, J1, K1, D3*D1*D2W, pm);
        }
    }
}

static double
cic_readout_tuned(FastPMPainter * painter, FastPMFloat * canvas, double pos[3])
{
    PM * pm = painter->pm;
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

    int I1=I+1; 
    int J1=J+1; 
    int K1=K+1; 

    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes

    while(UNLIKELY(I < 0)) I += pm->Nmesh[0];
    while(UNLIKELY(J < 0)) J += pm->Nmesh[1];
    while(UNLIKELY(K < 0)) K += pm->Nmesh[2];
    while(UNLIKELY(I >= pm->Nmesh[0])) I -= pm->Nmesh[0];
    while(UNLIKELY(J >= pm->Nmesh[1])) J -= pm->Nmesh[1];
    while(UNLIKELY(K >= pm->Nmesh[2])) K -= pm->Nmesh[2];


    while(UNLIKELY(I1 < 0)) I1 += pm->Nmesh[0];
    while(UNLIKELY(J1 < 0)) J1 += pm->Nmesh[1];
    while(UNLIKELY(K1 < 0)) K1 += pm->Nmesh[2];
    while(UNLIKELY(I1 >= pm->Nmesh[0])) I1-=pm->Nmesh[0];
    while(UNLIKELY(J1 >= pm->Nmesh[1])) J1-=pm->Nmesh[1];
    while(UNLIKELY(K1 >= pm->Nmesh[2])) K1-=pm->Nmesh[2];

    I -= pm->IRegion.start[0];
    J -= pm->IRegion.start[1];
    I1 -= pm->IRegion.start[0];
    J1 -= pm->IRegion.start[1];

    double value = 0;

    if(LIKELY(0 <= I && I < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                value += REd(canvas, I, J,  K,  T3*T1*T2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                value += REd(canvas, I, J,  K1, D3*T1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                value += REd(canvas, I, J1, K,  T3*T1*D2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                value += REd(canvas, I, J1, K1, D3*T1*D2W, pm);
        }
    }

    if(LIKELY(0 <= I1 && I1 < pm->IRegion.size[0])) {
        if(LIKELY(0 <= J && J < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                value += REd(canvas, I1, J,  K,  T3*D1*T2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                value += REd(canvas, I1, J,  K1, D3*D1*T2W, pm);
        }
        if(LIKELY(0 <= J1 && J1 < pm->IRegion.size[1])) {
            if(LIKELY(0 <= K && K < pm->IRegion.size[2]))
                value += REd(canvas, I1, J1, K,  T3*D1*D2W, pm);
            if(LIKELY(0 <= K1 && K1 < pm->IRegion.size[2]))
                value += REd(canvas, I1, J1, K1, D3*D1*D2W, pm);
        }
    }
    return value;
}

