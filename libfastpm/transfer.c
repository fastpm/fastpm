#include <math.h>
#include <stdlib.h>
#include <fastpm/libfastpm.h>
#include "pmpfft.h"
#include "transfer.h"

void 
fastpm_apply_smoothing_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double sml) 
{

#pragma omp parallel 
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        int d;
        int i;
        double *kernel[3];
        for(d = 0; d < 3; d ++) {
            kernel[d] = malloc(sizeof(double) * pm->Nmesh[d]);
            for(i = 0; i < pm->Nmesh[d]; i ++) {
                double kk = kiter.fac[d][i].kk;
                kernel[d][i] = exp(- 0.5 * kk * sml * sml);
            }
        }
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            for(dir = 0; dir < 3; dir++) 
                smth *= kernel[dir][kiter.iabs[dir]];

            /* - i k[d] */
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
        for(d = 0; d < 3; d ++) {
            free(kernel[d]);
        }
    }
}

void 
fastpm_apply_diff_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k_finite = kiter.fac[dir][kiter.iabs[dir]].k_finite;
            /* - i k[d] */
            to[kiter.ind + 0] =   from[kiter.ind + 1] * (k_finite);
            to[kiter.ind + 1] = - from[kiter.ind + 0] * (k_finite);
        }
    }
}

