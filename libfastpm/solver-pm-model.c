#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "pmghosts.h"
#include "pm2lpt.h"
#include "solver-pm-internal.h"

static double
measure_large_scale_power(FastPMModel * model, PMStore * p);

void fastpm_model_init(FastPMModel * model, FastPM * fastpm, FastPMModelType type)
{
    model->type = type;
    PMInit pminit = fastpm->pm_2lpt->init;
    if (fastpm->nc <= 256) {
        model->factor = 1;
        /* For LCDM cosmology, smaller mesh introduces too much aliasing noise. not good! */
        pminit.Nmesh = fastpm->nc * 2;
    } else {
        model->factor = 4;
        pminit.Nmesh = fastpm->nc / 2;
    }

    model->pm = malloc(sizeof(PM));
    model->fastpm = fastpm;
    pm_init(model->pm, &pminit, &fastpm->p->iface, fastpm->comm);
}

void fastpm_model_destroy(FastPMModel * model)
{
    pm_store_destroy(model->psub);
    free(model->psub);
    pm_destroy(model->pm);
    free(model->pm);
}

void fastpm_model_build(FastPMModel * model, PMStore * p, double ainit)
{
    PMStore * psub = malloc(sizeof(PMStore));
    pm_store_init(psub);
    pm_store_create_subsample(psub, p, PACK_POS | PACK_DX1 | PACK_DX2, model->factor, model->fastpm->nc);
    pm_2lpt_evolve(ainit, psub, model->fastpm->omega_m, 0);
    pm_store_wrap(psub, model->pm->BoxSize);
    model->Pexpect = measure_large_scale_power(model, psub);
    /* We no longer need to store the subsample once the initial linear power is calculated */
    model->a_x = ainit;
    model->psub = psub;
}

void fastpm_model_evolve(FastPMModel * model, double af)
{
    switch(model->type) {
        case FASTPM_MODEL_LINEAR:
        {
            model->Pexpect /= pow(fastpm_growth_factor(model->fastpm, model->a_x), 2);
            model->Pexpect *= pow(fastpm_growth_factor(model->fastpm, af), 2);
        }
        break;
        case FASTPM_MODEL_2LPT:
        {
            FastPMDrift drift;
            fastpm_drift_init(&drift, model->fastpm, model->psub, af);
            int i;
            for(i = 0; i < model->psub->np; i ++) {
                double xo[3];
                fastpm_drift_one_2lpt(&drift, i, xo);
                int d;
                for(d = 0; d < 3; d ++) {
                    model->psub->x[i][d] = xo[d];
                }
            }
            pm_store_wrap(model->psub, model->pm->BoxSize);
            model->psub->a_x = af;
            model->psub->a_v = af;
            model->Pexpect = measure_large_scale_power(model, model->psub);
        }
        break;

    }
    fastpm_info("Pexpect = %g, af=%g\n", model->Pexpect, af);
    model->a_x = af;
}

static void
scale_acc(PMStore * po, double correction, double fudge)
{
    ptrdiff_t i;
    correction = pow(correction, fudge);

#pragma omp parallel for
    for(i = 0; i < po->np; i ++) {
        int d;
        for(d = 0; d < 3; d ++) {
            po->acc[i][d] *= correction;
        }
    }
}


static double 
find_correction_eval(double correction, void * data)
{
    FastPMModel * model = (FastPMModel*) data;
    FastPM * fastpm = model->fastpm;
    PMStore * po = model->ev.po;

    scale_acc(po, correction, 1.0);

    fastpm_kick_store(fastpm, po, po, model->ev.a_v1);
    fastpm_drift_store(fastpm, po, po, model->ev.a_x1);

    double Plarge = measure_large_scale_power(model, po);

    fastpm_drift_store(fastpm, po, po, model->ev.a_x);
    fastpm_kick_store(fastpm, po, po, model->ev.a_v);
    scale_acc(po, 1.0 / correction, 1.0);

    double res = Plarge / model->Pexpect;
    return res - 1.0;
}

double
fastpm_model_find_correction(FastPMModel * model,
        double a_x, double a_x1, double a_v, double a_v1)
{

    if(a_x == a_x1) {
        /* no correction is needed in this case*/
        /* FIXME: when step size is sufficiently small, use 1.0 too */
        return 1.0;
    }
    PMStore * po = alloca(sizeof(PMStore));
    pm_store_init(po);

    /* Watch out: the same subsample as fastpm->pmodel */
    pm_store_create_subsample(po, model->fastpm->p,
            PACK_POS| PACK_VEL | PACK_ACC | PACK_DX1 | PACK_DX2, model->factor, model->fastpm->nc);

    gsl_root_fsolver * s;
    gsl_function F;

    F.function = find_correction_eval;
    F.params = (void*) model;
    model->ev.po = po;
    model->ev.a_x = a_x;
    model->ev.a_x1 = a_x1;
    model->ev.a_v = a_v;
    model->ev.a_v1 = a_v1;

    int iter = 0;
    double r_hi = 0, r_lo = 0;
    double x = 0;
    double x_lo = 0.9;
    double x_hi = 1.1;
    while((r_hi = find_correction_eval(x_hi, model)) < 0) {
        if(-r_hi < 1e-7) break;
        iter ++;
        fastpm_info("iter = %d x_hi = %g r = %g\n", iter, x_hi, r_hi);
        x_hi *= 1.2;
    }
    while((r_lo = find_correction_eval(x_lo, model)) > 0) {
        if(r_lo < 1e-7) break;
        iter ++;
        fastpm_info("iter = %d x_lo = %g r= %g\n", iter, x_lo, r_lo);
        x_lo /= 1.2;
    }

    if(fabs(r_hi) < 1e-7 && fabs(r_lo) < 1e-7) {
        pm_store_destroy(po);
        return 1.0;
    }

    s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    int status;
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        x = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                0, 1e-5);
        double r = find_correction_eval(x, model);
        fastpm_info("iter = %d x = %g r = %g\n", iter, x, r);
    }
    while (status == GSL_CONTINUE && iter < 10);
    gsl_root_fsolver_free(s);

    pm_store_destroy(po);

    return x;
}

static double 
measure_large_scale_power(FastPMModel * model, PMStore * p)
{
    PM * pm = model->pm;

    FastPMFloat * canvas = pm_alloc(pm);
    FastPMFloat * delta_k = pm_alloc(pm);

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    /* Note that pm_calculate_linear_power will divide by the 0-th mode
     * thus we do not need to set the mass of particles correctly */
    pm_paint(pm, canvas, p, p->np + pgd->nghosts, 1.0);

    pm_r2c(pm, canvas, delta_k);

    pm_ghosts_free(pgd);

    double Plin = pm_calculate_linear_power(pm, delta_k, model->fastpm->K_LINEAR);

    pm_free(pm, delta_k);
    pm_free(pm, canvas);
    return Plin;
}
