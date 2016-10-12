#include <string.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/transfer.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmghosts.h"

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
apply_force_kernel_eastwood(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir, int cic) 
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
    if(cic) {
        /* now sharpen for mass assignment */
        /* L1 */
        fastpm_apply_decic_transfer(pm, to, to);
        /* L2 */
        fastpm_apply_decic_transfer(pm, to, to);
    }
}

static void 
apply_force_kernel_gadget(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir, int cic)
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
    if(cic) {
        /* now sharpen for mass assignment */
        /* L1 */
        fastpm_apply_decic_transfer(pm, to, to);
        /* L2 */
        fastpm_apply_decic_transfer(pm, to, to);
    }
}

static void
apply_gaussian_dealiasing(PM * pm, FastPMFloat * from, FastPMFloat * to, double N)
{
    /* N is rms in mesh size */
    double r0 = N * pm->BoxSize[0] / pm->Nmesh[0];

#pragma omp parallel 
    {
        PMKIter kiter;
        int d;
        int i;
        double *kernel[3];
        pm_kiter_init(pm, &kiter);
        for(d = 0; d < 3; d ++) {
            kernel[d] = malloc(sizeof(double) * pm->Nmesh[d]);
            for(i = 0; i < pm->Nmesh[d]; i ++) {
                kernel[d][i] = exp(- 0.5 * pow(kiter.k[d][i] * r0, 2));
            }
        }

        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double fac = 1;
            ptrdiff_t ind = kiter.ind;
            for(d = 0; d < 3; d++) {
                fac *= kernel[d][kiter.iabs[d]];
            }
            to[ind + 0] *= fac;
            to[ind + 1] *= fac;
        }
    }
}
static double
gaussian36(double k, double * knq)
{
    double x = k / *knq;
    return exp(- 36 * pow(x, 36));
}

void
fastpm_gravity_calculate(FastPMGravity * gravity,
    PM * pm,
    FastPMStore * p,
    FastPMFloat * delta_k)
{
    FastPMPainter reader[1];
    FastPMPainter painter[1];

    fastpm_painter_init(reader, pm, FASTPM_PAINTER_CIC, 1);
    fastpm_painter_init(painter, pm, gravity->PainterType, gravity->PainterSupport);

    /* watch out: boost the density since mesh is finer than grid */
    long long np = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, pm_comm(pm));

    double density_factor = pm->Norm / np;

    CLOCK(ghosts);
    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);
    LEAVE(ghosts);

    FastPMFloat * canvas = pm_alloc(pm);

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * We thus have to boost the density by density_factor.
     *
     * This gives us over density + 1
     *
     * because rhobar = N_g ^3 / V
     * we paint rho V / (B N_g^3) * B = rho / rhobar. The last B is the extra density factor.
     * */
    CLOCK(paint);
    fastpm_paint_local(painter, canvas,
                p, p->np + pgd->nghosts, NULL, 0);
    fastpm_apply_multiply_transfer(pm, canvas, canvas, density_factor);
    LEAVE(paint);

    CLOCK(r2c);
    pm_r2c(pm, canvas, delta_k);
    LEAVE(r2c);

    /* calculate the forces save them to p->acc */

    switch(gravity->DealiasingType) {
        case FASTPM_DEALIASING_TWO_THIRD:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_lowpass_transfer(pm, delta_k, delta_k, 2.0 / 3 * k_nq);
            }
        break;
        case FASTPM_DEALIASING_GAUSSIAN:
            apply_gaussian_dealiasing(pm, delta_k, delta_k, 1.0);
        break;
        case FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN:
            apply_gaussian_dealiasing(pm, delta_k, delta_k, 4.0);
        break;
        case FASTPM_DEALIASING_GAUSSIAN36:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) gaussian36, &k_nq);
        }
        break;
        case FASTPM_DEALIASING_NONE:
        break;
        default:
            fastpm_raise(-1, "wrong dealiasing kernel type");
    }

    int d;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        CLOCK(transfer);
        switch(gravity->KernelType) {
            case FASTPM_KERNEL_EASTWOOD:
                apply_force_kernel_eastwood(pm, delta_k, canvas, d, 1);
            break;
            case FASTPM_KERNEL_NAIVE:
                apply_force_kernel_eastwood(pm, delta_k, canvas, d, 0);
            break;
            case FASTPM_KERNEL_GADGET:
                apply_force_kernel_gadget(pm, delta_k, canvas, d, 1);
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
        fastpm_readout_local(reader, canvas, p, p->np + pgd->nghosts, NULL, ACC[d]);
        LEAVE(readout);

        CLOCK(reduce);
        pm_ghosts_reduce(pgd, ACC[d]);
        LEAVE(reduce);
    }
    
    int i;
    
    for(i = 0; i < np; i++) {
        p->potential[i] = 0.0;
    }

    pm_free(pm, canvas);

    pm_ghosts_free(pgd);
}

