#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include <fastpm/libfastpm.h>
#include <fastpm/cosmology.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "pmghosts.h"
#include "pm2lpt.h"

static void 
apply_za_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{
    /* This is the force in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.fac[dir][kiter.iabs[dir]].k;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.fac[d][kiter.iabs[d]].kk;
            }
            ptrdiff_t ind = kiter.ind;
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] = - from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = + from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}


static void 
apply_2lpt_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir1, int dir2) 
{
    /* This is the force in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {

            int d;
            double k1_finite = kiter.fac[dir1][kiter.iabs[dir1]].k;
            double k2_finite = kiter.fac[dir2][kiter.iabs[dir2]].k;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.fac[d][kiter.iabs[d]].kk;
            }
            ptrdiff_t ind = kiter.ind;
            if(kk_finite > 0) {
                to[ind + 0] = from[ind + 0] * (-k1_finite * k2_finite / kk_finite);
                to[ind + 1] = from[ind + 1] * (-k1_finite * k2_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}


void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, PMStore * p, double shift[3]) 
{
/* calculate dx1, dx2, for initial fluctuation delta_k.
 * shift: martin has shift = 0.5, 0.5, 0.5.
 * Use shift of 0, 0, 0 if in doublt. 
 *   */

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    ptrdiff_t i;
    int d;

    FastPMFloat * workspace = pm_alloc(pm);
    FastPMFloat * source =  pm_alloc(pm);
    memset(source, 0, sizeof(source[0]) * pm->allocsize);

    FastPMFloat * field[3];

    for(d = 0; d < 3; d++ )
        field[d] = pm_alloc(pm);

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] -= shift[d];
        }
    }
     
    int DX1[] = {PACK_DX1_X, PACK_DX1_Y, PACK_DX1_Z};
    int DX2[] = {PACK_DX2_X, PACK_DX2_Y, PACK_DX2_Z};
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};

    for(d = 0; d < 3; d++) {

        apply_za_transfer(pm, delta_k, workspace, d);

        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {        
            p->dx1[i][d] = pm_readout_one(pm, workspace, p, i);
        }
        pm_ghosts_reduce(pgd, DX1[d]);
    } 

    for(d = 0; d< 3; d++) {
        apply_2lpt_transfer(pm, delta_k, field[d], d, d);
        pm_c2r(pm, field[d]);
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }    
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];

        apply_2lpt_transfer(pm, delta_k, workspace, d1, d2);
        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] -= workspace[i] * workspace[i];
        }
    } 
    pm_r2c(pm, source, workspace);
    pm_assign(pm, workspace, source);

    for(d = 0; d < 3; d++) {
        /* 
         * We absorb some the negative factor in za transfer to below;
         *
         * */
        apply_za_transfer(pm, source, workspace, d);
        pm_c2r(pm, workspace);

#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {        
            /* this ensures x = x0 + dx1(t) + 3/ 7 dx2(t) */
            p->dx2[i][d] = pm_readout_one(pm, workspace, p, i) / pm->Norm ;
        }
        pm_ghosts_reduce(pgd, DX2[d]);
    }

#ifdef PM_2LPT_DUMP
    fwrite(p->dx1, sizeof(p->dx1[0]), p->np, fopen("dx1.f4x3", "w"));
    fwrite(p->dx2, sizeof(p->dx2[0]), p->np, fopen("dx2.f4x3", "w"));
#endif

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += shift[d];
        }
    }

    for(d = 0; d < 3; d ++) {
        pm_free(pm, field[2-d]);
    }
    pm_free(pm, source);
    pm_free(pm, workspace);

    pm_ghosts_free(pgd);

}

// Interpolate position and velocity for snapshot at a=aout
void 
pm_2lpt_evolve(double aout, PMStore * p, double Omega)
{
    int np = p->np;

    Cosmology c = {
            .OmegaM = Omega,
            .OmegaLambda = 1 - Omega,
        };

    const float Dplus = GrowthFactor(aout, c);

    const double omega=OmegaA(aout, c);
    const double D2 = Dplus*Dplus*pow(omega/Omega, -1.0/143.0);
    const double D20 = pow(Omega, -1.0/143.0);
    

    float Dv=DprimeQ(aout, 1.0, c); // dD_{za}/dy
    float Dv2=GrowthFactor2v(aout, c);   // dD_{2lpt}/dy

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* Use the more accurate 2LPT dx2 term */
            p->dx2[i][d] *= 3.0 / 7.0 * D20;

            p->x[i][d] += Dplus * p->dx1[i][d] + D2 * p->dx2[i][d];

            p->v[i][d] = (p->dx1[i][d]*Dv + p->dx2[i][d]*Dv2);
        }
    }
    p->a_x = p->a_v = aout;
}

