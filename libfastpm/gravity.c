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
    fastpm_apply_laplace_transfer(pm, from, to, order);
    fastpm_apply_multiply_transfer(pm, to, to, -1);
}

static void
apply_grad_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir, int order)
{
    /* no need to print these, since we will check for FFT fields with pm_check_values.*/
    if(0) {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        float ** klist[2] = {kiter.k, kiter.k_finite};
        int i;
        for(i = 0; i < pm_nmesh(pm)[dir]; i ++) {
            double k_finite = klist[order][dir][i];
            fastpm_info("fourier space kernel[%d] = %g\n", i, k_finite);
        }
    }
#pragma omp parallel 
    {
        PMKIter kiter;
        ptrdiff_t * Nmesh = pm_nmesh(pm);
        float ** klist[2] = {kiter.k, kiter.k_finite};
        pm_kiter_init(pm, &kiter);
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
apply_gaussian_softening(PM * pm, FastPMFloat * from, FastPMFloat * to, double N)
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
fastpm_kernel_type_get_orders(FastPMKernelType type,
    int *potorder,
    int *gradorder,
    int *deconvolveorder)
{
    switch(type) {
        case FASTPM_KERNEL_EASTWOOD:
            *potorder = 0;
            *gradorder = 0;
            /* now sharpen for mass assignment */
            /* L1 and L2*/
            *deconvolveorder = 2;
        break;
        case FASTPM_KERNEL_NAIVE:
            *potorder = 0;
            *gradorder = 0;
            *deconvolveorder = 0;
        break;
        case FASTPM_KERNEL_GADGET:
            *potorder = 0;
            *gradorder = 1;
            *deconvolveorder = 2;
        break;
        case FASTPM_KERNEL_1_4:
            *potorder = 0;
            *gradorder = 1;
            *deconvolveorder = 0;
        break;
        case FASTPM_KERNEL_3_4:
            *potorder = 1;
            *gradorder = 1;
            *deconvolveorder = 0;
        break;
        case FASTPM_KERNEL_5_4:
            *potorder = 2;
            *gradorder = 1;
            *deconvolveorder = 0;
        break;
        case FASTPM_KERNEL_3_2:
            *potorder = 1;
            *gradorder = 0;
            *deconvolveorder = 0;
        break;
        default:
            fastpm_raise(-1, "Wrong kernel type\n");
    }
}

void
gravity_apply_kernel_transfer(FastPMKernelType type,
        PM * pm,
        FastPMFloat * delta_k,
        FastPMFloat * canvas, FastPMFieldDescr field)
{
    int potorder, gradorder, deconvolveorder;
    fastpm_kernel_type_get_orders(type, &potorder, &gradorder, &deconvolveorder);

    while(deconvolveorder > 0) {
        fastpm_apply_decic_transfer(pm, canvas, canvas);
        deconvolveorder--;
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
apply_softening_transfer(FastPMSofteningType type, PM * pm, FastPMFloat * from, FastPMFloat * to)
{
    switch(type) {
        case FASTPM_SOFTENING_TWO_THIRD:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_lowpass_transfer(pm, from, to, 2.0 / 3 * k_nq);
            }
        break;
        case FASTPM_SOFTENING_GAUSSIAN:
            apply_gaussian_softening(pm, from, to, 1.0);
        break;
        case FASTPM_SOFTENING_GADGET_LONG_RANGE:
            apply_gaussian_softening(pm, from, to, pow(2,0.5)*1.25);
        break;
        case FASTPM_SOFTENING_GAUSSIAN36:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_any_transfer(pm, from, to, (fastpm_fkfunc) gaussian36, &k_nq);
            }
        break;
        case FASTPM_SOFTENING_NONE:
        break;
        default:
            fastpm_raise(-1, "wrong softening kernel type");
    }
}

void
_fastpm_solver_create_ghosts(FastPMSolver * fastpm, int support, PMGhostData * pgd[6])
{
    PM * pm = fastpm->pm;

    CLOCK(ghosts);

    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;
        pgd[si] = pm_ghosts_create(pm, p, p->attributes, support);
        pm_ghosts_send(pgd[si], COLUMN_POS);
        pm_ghosts_send(pgd[si], COLUMN_ID);
        if(p->mass)
            pm_ghosts_send(pgd[si], COLUMN_MASS);
    }
    LEAVE(ghosts);
}

void
_fastpm_solver_destroy_ghosts(FastPMSolver * fastpm, PMGhostData * pgd[6])
{
    int si;
    for(si = FASTPM_SOLVER_NSPECIES - 1; si >= 0; si--) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;
        pm_ghosts_free(pgd[si]);
    }

}


void
_fastpm_solver_compute_delta_k(FastPMSolver * fastpm, FastPMPainter * painter, PMGhostData * pgd[6], FastPMFloat * canvas, FastPMFloat * delta_k)
{
    PM * pm = fastpm->pm;

    double total_mass = 0;

    FastPMFieldDescr FASTPM_FIELD_DESCR_NONE = {0, 0};

    pm_clear(pm, canvas);
    /* Watch out: paint paints the mass per cell;
     * divide by mean mass per cell to convert to matter overdensity, which
     * goes into Poisson's equation. 
     *
     * In this perspective, we are still operating with the dimension-less
     * Poisson's equation where the critical density factors canceled
     * with gravity constants into 1.5 OmegaM,
     * */
    CLOCK(paint);

    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;

        VALGRIND_CHECK_MEM_IS_DEFINED(p->x, sizeof(p->x[0]) * p->np);
        VALGRIND_CHECK_MEM_IS_DEFINED(pgd[si]->p->x, sizeof(pgd[si]->p->x[0]) * pgd[si]->p->np);

        double total_mass1 = 0;
        ptrdiff_t i;
        for (i = 0; i < p->np; i ++){
            total_mass1 += fastpm_store_get_mass(p, i);
        }
        total_mass += total_mass1;
        fastpm_paint_local(painter, canvas, p, p->np, FASTPM_FIELD_DESCR_NONE);
        fastpm_paint_local(painter, canvas, pgd[si]->p, pgd[si]->p->np, FASTPM_FIELD_DESCR_NONE);
    }
    LEAVE(paint);

    MPI_Allreduce(MPI_IN_PLACE, &total_mass, 1, MPI_DOUBLE, MPI_SUM, fastpm->comm);
    double mean_mass_per_cell = total_mass / pm->Norm;

    CLOCK(transfer);
    fastpm_apply_multiply_transfer(pm, canvas, canvas, 1.0 / mean_mass_per_cell);
    LEAVE(transfer);

    CLOCK(r2c);

    pm_check_values(pm, canvas, "After painting");
    pm_r2c(pm, canvas, delta_k);
    pm_check_values(pm, delta_k, "After r2c");

    LEAVE(r2c);

}

void
_fastpm_solver_compute_force(FastPMSolver * fastpm,
    FastPMPainter * reader,
    FastPMKernelType kernel,
    PMGhostData * pgd[6],
    FastPMFloat * canvas,
    FastPMFloat * delta_k, FastPMFieldDescr * ACC, int nacc)
{
    int d;
    PM * pm = fastpm->pm;

    CLOCK(transfer);
    CLOCK(c2r);
    CLOCK(readout);
    CLOCK(reduce);

    for(d = 0; d < nacc; d ++) {

        ENTER(transfer);
        gravity_apply_kernel_transfer(kernel, pm, delta_k, canvas, ACC[d]);
        LEAVE(transfer);

        ENTER(c2r);
        pm_check_values(pm, canvas, "Before c2r %d", d);
        pm_c2r(pm, canvas);
        pm_check_values(pm, canvas, "After c2r %d", d);
        LEAVE(c2r);

        ENTER(readout);
        int si;
        for(si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
            FastPMStore * p = fastpm_solver_get_species(fastpm, si);
            if(!p) continue;
            fastpm_readout_local(reader, canvas, p, p->np, ACC[d]);
            fastpm_readout_local(reader, canvas, pgd[si]->p, pgd[si]->p->np, ACC[d]);
        }
        LEAVE(readout);

    }

    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;
        double acc_std[3], acc_mean[3], acc_min[3], acc_max[3];
        fastpm_store_summary(p, COLUMN_ACC, pm_comm(pm), "<s->", acc_min, acc_std, acc_mean, acc_max);
        for(d = 0; d < 3; d ++) {
            fastpm_info("p%s    acc[%d]: %g %g %g %g\n",
                p->name, d, acc_min[d], acc_std[d], acc_mean[d], acc_max[d]);
        }
        fastpm_store_summary(pgd[si]->p, COLUMN_ACC, pm_comm(pm), "<s->", acc_min, acc_std, acc_mean, acc_max);
        for(d = 0; d < 3; d ++) {
            fastpm_info("ghost acc[%d]: %g %g %g %g\n",
                d, acc_min[d], acc_std[d], acc_mean[d], acc_max[d]);
        }
        fastpm_store_summary(p, COLUMN_ACC, pm_comm(pm), "<s->", acc_min, acc_std, acc_mean, acc_max);
        for(d = 0; d < 3; d ++) {
            fastpm_info("p%s+g  acc[%d]: %g %g %g %g\n",
                p->name, d, acc_min[d], acc_std[d], acc_mean[d], acc_max[d]);
        }

        ENTER(reduce);
        pm_ghosts_reduce(pgd[si], COLUMN_ACC, FastPMReduceAddFloat, NULL);

        if(p->potential != NULL) {
            pm_ghosts_reduce(pgd[si], COLUMN_POTENTIAL, FastPMReduceAddFloat, NULL);
        }
        LEAVE(reduce);
    }


}

void
fastpm_solver_compute_force(FastPMSolver * fastpm,
    FastPMPainter * painter,
    FastPMSofteningType dealias,
    FastPMKernelType kernel,
    FastPMFloat * delta_k)
{
    PM * pm = fastpm->pm;
    PMGhostData * pgd[FASTPM_SOLVER_NSPECIES];

    FastPMFloat * canvas = pm_alloc(pm);

    _fastpm_solver_create_ghosts(fastpm, painter->support, pgd);

    _fastpm_solver_compute_delta_k(fastpm, painter, pgd, canvas, delta_k);

    CLOCK(dealias);
    /* calculate the forces save them to p->acc */
    apply_softening_transfer(dealias, pm, delta_k, delta_k);
    pm_check_values(pm, delta_k, "After softening");
    LEAVE(dealias);

    FastPMFieldDescr ACC[] = {
         {COLUMN_ACC, 0},
         {COLUMN_ACC, 1},
         {COLUMN_ACC, 2},
         {COLUMN_POTENTIAL, 0}
    };

    int nacc = 4;
    /* skip potential if not wanted */

    if(NULL == fastpm_solver_get_species(fastpm, FASTPM_SPECIES_CDM)->potential) {
        nacc = 3;
    }

    _fastpm_solver_compute_force(fastpm, painter, kernel, pgd, canvas, delta_k, ACC, nacc);

    _fastpm_solver_destroy_ghosts(fastpm, pgd);

    pm_free(pm, canvas);
}
