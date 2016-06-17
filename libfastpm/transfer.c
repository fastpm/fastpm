#include <math.h>
#include <stdlib.h>
#include <fastpm/libfastpm.h>
#include "pmpfft.h"

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
                double kk = kiter.kk[d][i];
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
fastpm_apply_lowpass_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double kth)
{
    double kth2 = kth * kth;
#pragma omp parallel 
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            double kk = 0;
            for(dir = 0; dir < 3; dir++) {
                kk += kiter.kk[dir][kiter.iabs[dir]];
            }
            if(kk < kth2) smth = 1;
            else smth = 0;
            /* - i k[d] */
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
    }
}

static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}


void 
fastpm_apply_decic_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to) 
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
                double w = kiter.k[d][i] * pm->BoxSize[d] / pm->Nmesh[d];
                double cic = sinc_unnormed(0.5 * w);
                /* Watchout: this does divide by sinc, not sinc 2, */
                kernel[d][i] = 1.0 / pow(cic, 2);
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
            double k_finite = kiter.k_finite[dir][kiter.iabs[dir]];
            /* i k[d] */
            FastPMFloat tmp[2];
            tmp[0] = - from[kiter.ind + 1] * (k_finite);
            tmp[1] =   from[kiter.ind + 0] * (k_finite);
            to[kiter.ind + 0] = tmp[0];
            to[kiter.ind + 1] = tmp[1];
        }
    }
}

void fastpm_apply_laplace_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to) {
#pragma omp parallel
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk_finite = 0.;
            for(d = 0; d < 3; d++) {
                /* Referee says 1 / kk is better than 1 / sin(kk);
                 *
                 * In reality kk_finite seems to give a better agreement on
                 * HMC derivatives. The difference is ~ 0.1%
                 *
                 * On large scales this does not matter.
                 * */
                /* 2-point Finite differentiation */
                kk_finite += kiter.kk_finite[d][kiter.iabs[d]];
            }
            if(kk_finite == 0)
            {
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            }
            else
            {
                /* 1 / k**2 */
                to[kiter.ind + 0] =  from[kiter.ind + 0]  / kk_finite;
                to[kiter.ind + 1] =  from[kiter.ind + 1]  / kk_finite;
            }
        }
    }
}

void
fastpm_apply_any_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, fastpm_fkfunc func, void * data)
{
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            double kk = 0;
            for(dir = 0; dir < 3; dir++) {
                kk += kiter.kk[dir][kiter.iabs[dir]];
            }
            double k = sqrt(kk);
            smth = func(k, data);
            /* - i k[d] */
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
    }
}

void
fastpm_apply_multiply_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double value)
{
    ptrdiff_t i;
#pragma omp parallel for
    for(i = 0; i < pm_allocsize(pm); i ++) {
        to[i] = from[i] * value;
    }
}
