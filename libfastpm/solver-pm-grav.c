#include <string.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/transfer.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmghosts.h"
#include "pmstore.h"

static void 
apply_force_kernel(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{
    /* This is the potential gradient in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;

        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.k_finite[dir][kiter.iabs[dir]];
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.kk_finite[d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] =   from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = - from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

static void 
apply_force_kernel_5_4(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{
    /* This is the potential gradient in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;

        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.k_finite[dir][kiter.iabs[dir]];
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.kk_finite2[d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] =   from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = - from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

static void 
apply_force_kernel_3_2(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{
    /* This is the potential gradient in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;

        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.k[dir][kiter.iabs[dir]];
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.kk_finite[d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] =   from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = - from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

static void 
apply_force_kernel_eastwood(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
{
    /* This is the potential gradient in fourier space. - i k[dir] / k2 */

#pragma omp parallel 
    {
        PMKIter kiter;

        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.k[dir][kiter.iabs[dir]];
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.kk[d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] =   from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = - from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
    /* now sharpen for mass assignment */
    /* L1 */
    fastpm_apply_decic_transfer(pm, to, to);
    /* L2 */
    fastpm_apply_decic_transfer(pm, to, to);
}

void
fastpm_calculate_forces(FastPM * fastpm, FastPMFloat * delta_k)
{
    PMStore * p = fastpm->p;
    PM * pm = fastpm->pm;

    /* watch out: boost the density since mesh is finer than grid */
    double density_factor = fastpm->pm->Norm / pow(1.0 * fastpm->nc, 3);

    CLOCK(ghosts);
    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL); 
    LEAVE(ghosts);

    FastPMFloat * canvas = pm_alloc(pm);

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * We thus have to boost the density by density_factor.
     * */
    CLOCK(paint);
    pm_paint(pm, canvas, p, p->np + pgd->nghosts, density_factor);
    LEAVE(paint);

    CLOCK(r2c);
    pm_r2c(pm, canvas, delta_k);
    LEAVE(r2c);

    /* calculate the forces save them to p->acc */

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        CLOCK(transfer);
        switch(fastpm->KERNEL_TYPE) {
            case FASTPM_KERNEL_EASTWOOD:
                apply_force_kernel_eastwood(pm, delta_k, canvas, d);
            break;
            case FASTPM_KERNEL_3_4:
                apply_force_kernel(pm, delta_k, canvas, d);
            break;
            case FASTPM_KERNEL_5_4:
                apply_force_kernel_5_4(pm, delta_k, canvas, d);
            break;
            case FASTPM_KERNEL_3_2:
                apply_force_kernel_3_2(pm, delta_k, canvas, d);
            break;
            default:
                fastpm_raise(-1, "Wrong kernel type\n");
        }
        LEAVE(transfer);

        CLOCK(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        CLOCK(readout);
#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {
            p->acc[i][d] = pm_readout_one(pm, canvas, p, i);
        }
        LEAVE(readout);

        CLOCK(reduce);
        pm_ghosts_reduce(pgd, ACC[d]);
        LEAVE(reduce);
    }

    pm_free(pm, canvas);

    pm_ghosts_free(pgd);
}

