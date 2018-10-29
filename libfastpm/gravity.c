#include <string.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/transfer.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmghosts.h"

static void
apply_pot_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int order)
{
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        float ** kklist [3] = {kiter.kk, kiter.kk_finite, kiter.kk_finite2};
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kklist[order][d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - 1 / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] = - from[ind + 0] * (1 / kk_finite);
                to[ind + 1] = - from[ind + 1] * (1 / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

static void
apply_grad_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir, int order)
{
#pragma omp parallel 
    {
        PMKIter kiter;
        ptrdiff_t * Nmesh = pm_nmesh(pm);
        pm_kiter_init(pm, &kiter);
        float ** klist[2] = {kiter.k, kiter.k_finite};
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k_finite;
            k_finite = klist[order][dir][kiter.iabs[dir]];
            ptrdiff_t ind = kiter.ind;
            /* i k[d] */
            /* Watch out the data dependency */
            if(
                kiter.iabs[0] == (Nmesh[0] - kiter.iabs[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - kiter.iabs[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - kiter.iabs[2]) % Nmesh[2]
            ) {
                /* We are at the nyquist and the diff operator shall be zero;
                 * otherwise the force is not real! */
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            } else {
                FastPMFloat tmp = from[ind + 0] * (k_finite);
                to[ind + 0] = - from[ind + 1] * (k_finite);
                to[ind + 1] = tmp;
            }
        }
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
        for(d = 0; d < 3; d ++) {
            free(kernel[d]);
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
gravity_apply_kernel_transfer(FastPMGravity * gravity,
        PM * pm,
        FastPMFloat * delta_k,
        FastPMFloat * canvas, FastPMFieldDescr field)
{
    int potorder = 0;
    int gradorder = 0;
    switch(gravity->KernelType) {
        case FASTPM_KERNEL_EASTWOOD:
            potorder = 0;
            gradorder = 0;
            /* now sharpen for mass assignment */
            /* L1 */
            fastpm_apply_decic_transfer(pm, canvas, canvas);
            /* L2 */
            fastpm_apply_decic_transfer(pm, canvas, canvas);
        break;
        case FASTPM_KERNEL_NAIVE:
            potorder = 0;
            gradorder = 0;
        break;
        case FASTPM_KERNEL_GADGET:
            potorder = 0;
            gradorder = 1;
        break;
        case FASTPM_KERNEL_3_4:
            potorder = 1;
            gradorder = 1;
        break;
        case FASTPM_KERNEL_5_4:
            potorder = 2;
            gradorder = 1;
        break;
        case FASTPM_KERNEL_3_2:
            potorder = 1;
            gradorder = 0;
        break;
        default:
            fastpm_raise(-1, "Wrong kernel type\n");
    }

    int d1, d2;

    switch(field.attribute) {
        case COLUMN_POTENTIAL:
            apply_pot_transfer(pm, delta_k, canvas, potorder);
            break;
        case COLUMN_DENSITY:
            fastpm_apply_multiply_transfer(pm, delta_k, canvas, 1.0);
            break;
        case COLUMN_TIDAL:
            switch(field.memb) {
                case 0:
                    d1 = 0; d2 = 0;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
                case 1:
                    d1 = 1; d2 = 1;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
                case 2:
                    d1 = 2; d2 = 2;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
                case 3:
                    d1 = 0; d2 = 1;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
                case 4:
                    d1 = 1; d2 = 2;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
                case 5:
                    d1 = 2; d2 = 0;
                    apply_pot_transfer(pm, delta_k, canvas, potorder);
                    apply_grad_transfer(pm, canvas, canvas, d1, gradorder);
                    apply_grad_transfer(pm, canvas, canvas, d2, gradorder);
                    break;
            }
            break;
        case COLUMN_ACC:
            apply_pot_transfer(pm, delta_k, canvas, potorder);
            apply_grad_transfer(pm, canvas, canvas, field.memb, gradorder);
            break;
        default:
            fastpm_raise(-1, "Unknown type for gravity attribute\n");
    }
}
static void
apply_dealiasing_transfer(FastPMGravity * gravity, PM * pm, FastPMFloat * from, FastPMFloat * to)
{
    switch(gravity->DealiasingType) {
        case FASTPM_DEALIASING_TWO_THIRD:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_lowpass_transfer(pm, from, to, 2.0 / 3 * k_nq);
            }
        break;
        case FASTPM_DEALIASING_GAUSSIAN:
            apply_gaussian_dealiasing(pm, from, to, 1.0);
        break;
        case FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN:
            apply_gaussian_dealiasing(pm, from, to, 4.0);
        break;
        case FASTPM_DEALIASING_GAUSSIAN36:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_any_transfer(pm, from, to, (fastpm_fkfunc) gaussian36, &k_nq);
            }
        break;
        case FASTPM_DEALIASING_NONE:
        break;
        default:
            fastpm_raise(-1, "wrong dealiasing kernel type");
    }
}

void
fastpm_gravity_calculate(FastPMGravity * gravity,
    PM * pm,
    FastPMStore * p,
    FastPMFloat * delta_k)
{
    FastPMPainter reader[1];
    FastPMPainter painter[1];

    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);
    fastpm_painter_init(painter, pm, gravity->PainterType, gravity->PainterSupport);

    /* watch out: boost the density since mesh is finer than grid */
    long long np = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, pm_comm(pm));

    double density_factor = pm->Norm / np;

    CLOCK(ghosts);
    PMGhostData * pgd = pm_ghosts_create(pm, p, p->attributes);
    pm_ghosts_send(pgd, COLUMN_POS);
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

    VALGRIND_CHECK_MEM_IS_DEFINED(p->x, sizeof(p->x[0]) * p->np);
    VALGRIND_CHECK_MEM_IS_DEFINED(pgd->p->x, sizeof(pgd->p->x[0]) * pgd->p->np);
    FastPMFieldDescr FASTPM_FIELD_DESCR_NONE = {0, 0};

    pm_clear(pm, canvas);
    fastpm_paint_local(painter, canvas, p, p->np, FASTPM_FIELD_DESCR_NONE);
    fastpm_paint_local(painter, canvas, pgd->p, pgd->p->np, FASTPM_FIELD_DESCR_NONE);

    fastpm_apply_multiply_transfer(pm, canvas, canvas, density_factor);

    LEAVE(paint);
    CLOCK(r2c);
    pm_r2c(pm, canvas, delta_k);
    LEAVE(r2c);

    /* calculate the forces save them to p->acc */
    apply_dealiasing_transfer(gravity, pm, delta_k, delta_k);

    int d;

    FastPMFieldDescr ACC[] = {
         {COLUMN_ACC, 0},
         {COLUMN_ACC, 1},
         {COLUMN_ACC, 2},
         {COLUMN_POTENTIAL, 0}
    };

    CLOCK(transfer);
    CLOCK(c2r);
    CLOCK(readout);
    CLOCK(reduce);

    for(d = 0; d < 3 + 1; d ++) {
        /* skip potential if not wanted */
        if(p->potential == NULL && ACC[d].attribute == COLUMN_POTENTIAL) continue;
        ENTER(transfer);
        gravity_apply_kernel_transfer(gravity, pm, delta_k, canvas, ACC[d]);
        LEAVE(transfer);

        ENTER(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        ENTER(readout);
        fastpm_readout_local(reader, canvas, p, p->np, ACC[d]);
        fastpm_readout_local(reader, canvas, pgd->p, pgd->p->np, ACC[d]);
        LEAVE(readout);

    }

    ENTER(reduce);
    pm_ghosts_reduce(pgd, COLUMN_ACC, FastPMReduceAddFloat, NULL);

    if(p->potential != NULL) {
        pm_ghosts_reduce(pgd, COLUMN_POTENTIAL, FastPMReduceAddFloat, NULL);
    }

    LEAVE(reduce);

    pm_free(pm, canvas);

    pm_ghosts_free(pgd);
}
