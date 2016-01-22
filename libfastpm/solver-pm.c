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

#include "pmpfft.h"
#include "pmstore.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"

static void 
fastpm_set_time(FastPM * fastpm, 
    int istep,
    double * time_step,
    int nstep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1);

static void
fastpm_decompose(FastPM * fastpm);

void 
fastpm_kick(FastPM * fastpm, 
              PMStore * pi, PMStore * po,
              double af, int verbose);

void 
fastpm_drift(FastPM * fastpm,
               PMStore * pi, PMStore * po,
               double af, int verbose);

void 
fastpm_set_snapshot(FastPM * fastpm,
                PMStore * p, PMStore * po,
                double aout);

double fastpm_growth_factor(FastPM * fastpm, double a);


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void
scale_acc(PMStore * po, double correction, double fudge);
static double
measure_linear_power(FastPM * fastpm, PMStore * p, double a_x);
static double
find_correction(FastPM * fastpm, double Plin, 
        double a_x, double a_x1, double a_v, double a_v1);

/* Useful stuff */
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);

void fastpm_init(FastPM * fastpm, 
    int NprocY, 
    int UseFFTW, 
    MPI_Comm comm) {

    PMInit baseinit = {
            .Nmesh = fastpm->nc,
            .BoxSize = fastpm->boxsize, 
            .NprocY = NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = UseFFTW,
        };

    PMInit pm_2lptinit = baseinit;
    PMInit pm_linearinit = baseinit;

    pm_linearinit.Nmesh = fastpm->nc / 4;

    fastpm->comm = comm;
    MPI_Comm_rank(comm, &fastpm->ThisTask);
    MPI_Comm_size(comm, &fastpm->NTask);

    fastpm->p= malloc(sizeof(PMStore));
    fastpm->pm_2lpt = malloc(sizeof(PM));
    fastpm->pm_linear = malloc(sizeof(PM));

    pm_store_init(fastpm->p);

    pm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC, 
        fastpm->alloc_factor, comm);

    fastpm->vpm_list = vpm_create(fastpm->vpminit,
                           &baseinit, &fastpm->p->iface, comm);

    pm_init(fastpm->pm_2lpt, &pm_2lptinit, &fastpm->p->iface, comm);
    pm_init(fastpm->pm_linear, &pm_linearinit, &fastpm->p->iface, comm);

    int i = 0;
    for (i = 0; i < FASTPM_EXT_MAX; i ++) {
        fastpm->exts[i] = NULL;
    }

}

void 
fastpm_setup_ic(FastPM * fastpm, FastPMFloat * delta_k_ic, double ainit)
{
    if(delta_k_ic) {
        double shift[3] = {
            fastpm->boxsize / fastpm->nc * 0.5,
            fastpm->boxsize / fastpm->nc * 0.5,
            fastpm->boxsize / fastpm->nc * 0.5,
            };

        pm_store_set_lagrangian_position(fastpm->p, fastpm->pm_2lpt, shift);

        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(fastpm->pm_2lpt, delta_k_ic, fastpm->p, shift);
    }
    if(fastpm->USE_DX1_ONLY == 1) {
        memset(fastpm->p->dx2, 0, sizeof(fastpm->p->dx2[0]) * fastpm->p->np);
    }

    pm_store_summary(fastpm->p, fastpm->comm);
    pm_2lpt_evolve(ainit, fastpm->p, fastpm->omega_m, fastpm->USE_DX1_ONLY);
}

void
fastpm_evolve(FastPM * fastpm, double * time_step, int nstep) 
{
    CLOCK(decompose);
    CLOCK(force);
    CLOCK(kick);
    CLOCK(drift);
    CLOCK(afterforce);
    CLOCK(afterkick);
    CLOCK(afterdrift);
    CLOCK(correction);

    FastPMExtension * ext;

    MPI_Barrier(fastpm->comm);

    CLOCK(warmup);

    if(fastpm->USE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(fastpm->p->v, 0, sizeof(fastpm->p->v[0]) * fastpm->p->np);
    }

    LEAVE(warmup);

    MPI_Barrier(fastpm->comm);



    double correction = 1.0;
    double Plin0sub = 0;
    double Plin0 = 0;
    /* The last step is the 'terminal' step */
    int istep;
    for (istep = 0; istep < nstep; istep++) {
        double a_v, a_x, a_v1, a_x1;

        /* begining and ending of drift(x) and kick(v)*/
        fastpm_set_time(fastpm, istep, time_step, nstep,
                    &a_x, &a_x1, &a_v, &a_v1);

        fastpm_info("==== Step %d a_x = %6.4f a_x1 = %6.4f a_v = %6.4f a_v1 = %6.4f Nmesh = %d ====\n", 
                    istep, a_x, a_x1, a_v, a_v1, fastpm->pm->init.Nmesh);

        ENTER(decompose);
        fastpm_decompose(fastpm);
        LEAVE(decompose);

        /* Calculate PM forces. */
        FastPMFloat * delta_k = pm_alloc(fastpm->pm);

        /* watch out: boost the density since mesh is finer than grid */
        double density_factor = fastpm->pm->Norm / pow(1.0 * fastpm->nc, 3);

        ENTER(force);
        pm_calculate_forces(fastpm->p, fastpm->pm, delta_k, density_factor);
        LEAVE(force);

        double Plin = pm_calculate_linear_power(fastpm->pm, delta_k, fastpm->K_LINEAR);
        Plin /= pow(fastpm_growth_factor(fastpm, a_x), 2.0);

        if(istep == 0) {
            Plin0 = Plin;
            PMStore * po = alloca(sizeof(PMStore));
            pm_store_init(po);
            pm_store_create_subsample(po, fastpm->p, PACK_POS| PACK_VEL | PACK_ACC | PACK_ID, 4, fastpm->nc);
            Plin0sub = measure_linear_power(fastpm, po, a_x);
            pm_store_destroy(po);
        }

        fastpm_info("Last Step: <P(k<%g)> = %g Linear Theory = %g, correction=%g, res=%g\n", 
                          fastpm->K_LINEAR * 6.28 / fastpm->boxsize, Plin, Plin0, correction, Plin / Plin0); 

        ENTER(afterforce);
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_FORCE];
            ext; ext = ext->next) {
                ((fastpm_ext_after_force) ext->function) 
                    (fastpm, delta_k, a_x, ext->userdata);
        }
        LEAVE(afterforce);


        pm_free(fastpm->pm, delta_k);


        ENTER(afterdrift);
        /* take snapshots if needed, before the kick */
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_DRIFT];
            ext; ext = ext->next) {
                ((fastpm_ext_after_drift) ext->function) 
                    (fastpm, ext->userdata);
        }
        LEAVE(afterdrift);

        /* never go beyond 1.0 */
        if(a_x >= 1.0) break; 

        /* correct for linear theory before kick and drift */
        ENTER(correction);
        if(fastpm->USE_LINEAR_THEORY) {

            correction = find_correction(fastpm, Plin0sub, 
                        a_x, a_x1, a_v, a_v1);

            scale_acc(fastpm->p, correction, 1.0);
            /* the value of correction is leaked to the next step for reporting. */
        }
        LEAVE(correction);
        
        // Leap-frog "kick" -- velocities updated

        ENTER(kick);
        fastpm_kick(fastpm, fastpm->p, fastpm->p, a_v1, 1);
        LEAVE(kick);

        ENTER(afterkick);
        /* take snapshots if needed, before the drift */
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_KICK];
            ext; ext = ext->next) {
                ((fastpm_ext_after_kick) ext->function) 
                    (fastpm, ext->userdata);
        }
        LEAVE(afterkick);
        
        // Leap-frog "drift" -- positions updated

        ENTER(drift);
        fastpm_drift(fastpm, fastpm->p, fastpm->p, a_x1, 1);
        LEAVE(drift);


    }

}

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct find_correction_params {
    FastPM * fastpm;
    PMStore * po;
    double Plin0;
    double a_x;
    double a_x1;
    double a_v;
    double a_v1;
};

static double 
find_correction_eval(double correction, void * params) 
{
    struct find_correction_params * p = 
        (struct find_correction_params*) params;

    FastPM * fastpm = p->fastpm;
    PMStore * po = p->po;

    scale_acc(po, correction, 1.0);

    fastpm_kick(fastpm, po, po, p->a_v1, 0);
    fastpm_drift(fastpm, po, po, p->a_x1, 0);

    double Plin = measure_linear_power(fastpm, po, p->a_x1);
    fastpm_drift(fastpm, po, po, p->a_x, 0);
    fastpm_kick(fastpm, po, po, p->a_v, 0);
    scale_acc(po, 1.0 / correction, 1.0);

    double res = Plin / p->Plin0;
    return res - 1.0;
}

static double
find_correction(FastPM * fastpm, double Plin, 
        double a_x, double a_x1, double a_v, double a_v1) 
{

    PMStore * po = alloca(sizeof(PMStore));
    pm_store_init(po);
    pm_store_create_subsample(po, fastpm->p, PACK_POS| PACK_VEL | PACK_ACC | PACK_ID, 4, fastpm->nc);

    gsl_root_fsolver * s;
    gsl_function F;
    struct find_correction_params params = {
        fastpm, po, Plin, a_x, a_x1, a_v, a_v1
    };

    F.function = find_correction_eval;
    F.params = &params;

    s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    int iter = 0;
    double r = 0;
    double x_lo = 0.9;
    double x_hi = 1.1;

    while(find_correction_eval(x_hi, &params) < 0) {
        iter ++;
        fastpm_info("iter = %d x_hi = %g\n", iter, x_hi);
        x_hi *= 1.1;
        if(iter > 30) {
            fastpm_raise(-1, "Too many iterations finding lower bound\n");
        }
    }
    while(find_correction_eval(x_lo, &params) > 0) {
        iter ++;
        fastpm_info("iter = %d x_lo = %g\n", iter, x_lo);
        x_lo /= 1.1;
        if(iter > 60) {
            fastpm_raise(-1, "Too many iterations finding upper bound\n");
        }
    }
    int status;
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                0, 1e-2);
        fastpm_info("iter = %d correction = %g\n", iter, r);
    }
    while (status == GSL_CONTINUE && iter < 10);
    gsl_root_fsolver_free(s);

    pm_store_destroy(po);

    return r;
}

static double 
measure_linear_power(FastPM * fastpm, PMStore * p, double a_x) 
{
    PM * pm = fastpm->pm_linear;

    FastPMFloat * canvas = pm_alloc(pm);
    FastPMFloat * delta_k = pm_alloc(pm);

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    /* Note that pm_calculate_linear_power will divide by the 0-th mode
     * thus we do not need to set the mass of particles correctly */
    pm_paint(pm, canvas, p, p->np + pgd->nghosts, 1.0);

    pm_r2c(pm, canvas, delta_k);

    pm_ghosts_free(pgd);

    double Plin = pm_calculate_linear_power(pm, delta_k, fastpm->K_LINEAR);

    /* normalize to z=0 for comparison */
    Plin /= pow(fastpm_growth_factor(fastpm, a_x), 2.0);

    pm_free(pm, delta_k);
    pm_free(pm, canvas);
    return Plin;
}
void 
fastpm_destroy(FastPM * fastpm) 
{
    pm_store_destroy(fastpm->p);
    vpm_free(fastpm->vpm_list);
    pm_destroy(fastpm->pm_2lpt);
    pm_destroy(fastpm->pm_linear);
    free(fastpm->pm_linear);
    free(fastpm->pm_2lpt);
    /* FIXME: free VPM and stuff. */
    FastPMExtension * ext, * e2;
    int i;
    for(i = 0; i < FASTPM_EXT_MAX; i++) {
        for(ext = fastpm->exts[i]; ext; ext = e2) {
            e2 = ext->next;
            free(e2);
        }
    }
}

static int 
to_rank(void * pdata, ptrdiff_t i, void * data) 
{
    PMStore * p = (PMStore *) pdata;
    PM * pm = (PM*) data;
    double pos[3];
    p->iface.get_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
}

static void
fastpm_decompose(FastPM * fastpm) {
    /* apply periodic boundary and move particles to the correct rank */
    pm_store_wrap(fastpm->p, fastpm->pm->BoxSize);
    pm_store_decompose(fastpm->p, to_rank, fastpm->pm, fastpm->comm);
    size_t np_max;
    size_t np_min;
    /* FIXME move NTask to somewhere else. */
    double np_mean = pow(fastpm->nc, 3) / fastpm->pm->NTask;
    MPI_Allreduce(&fastpm->p->np, &np_max, 1, MPI_LONG, MPI_MAX, fastpm->comm);
    MPI_Allreduce(&fastpm->p->np, &np_min, 1, MPI_LONG, MPI_MIN, fastpm->comm);

    fastpm_info("Load imbalance is - %g / + %g\n",
        np_min / np_mean, np_max / np_mean);

}

#if 0
static void 
smooth_density(PM * pm, double r_s) 
{
    /* 
     *  This function smooth density by scale r_s. There could be a factor of sqrt(2)
     *  It is not used. */

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);
    {
        /* fill in the extra 'smoothing kernels' we will take the product */
        ptrdiff_t ind;
        int d;
        for(d = 0; d < 3; d++)
        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            fac[d][ind].extra = exp(-0.5 * fac[d][ind].kk * r_s * r_s);
        }
    }

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double smth = 1.;
            double kk = 0.;
            for(d = 0; d < 3; d++) {
                smth *= fac[d][i[d] + pm->ORegion.start[d]].extra;
                kk += fac[d][i[d] + pm->ORegion.start[d]].kk;
            }
            /* - i k[d] / k2 */
            if(LIKELY(kk> 0)) {
                pm->workspace[ind + 0] = pm->canvas[ind + 0] * smth;
                pm->workspace[ind + 1] = pm->canvas[ind + 1] * smth;
            } else {
                pm->workspace[ind + 0] = 0;
                pm->workspace[ind + 1] = 0;
            }
            pm_inc_o_index(pm, i);
        }
    }

    pm_destroy_k_factors(pm, fac);
}
#endif

void 
fastpm_interp(FastPM * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata) 
{
    /* interpolate and write snapshots, assuming p 
     * is at time a_x and a_v. */
    double a_x = fastpm->p->a_x;
    double a_v = fastpm->p->a_v;
    int iout;
    for(iout = 0; iout < nout; iout ++) {
        if(
        ! /* after a kick */
        (a_x < aout[iout] && aout[iout] < a_v)
        &&
        ! /* after a drift */
        (a_x >= aout[iout] && aout[iout] >= a_v)
        ) continue;

        PMStore * snapshot = alloca(sizeof(PMStore));
        pm_store_init(snapshot);
        pm_store_alloc(snapshot, fastpm->p->np_upper, PACK_ID | PACK_POS | PACK_VEL);

        fastpm_info("Taking a snapshot...\n");

        fastpm_set_snapshot(fastpm, fastpm->p, snapshot, aout[iout]);

        action(fastpm, snapshot, aout[iout], userdata);

        pm_store_destroy(snapshot);

    }
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

static void 
fastpm_set_time(FastPM * fastpm, 
    int istep,
    double * time_step,
    int nstep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1) 
{
    /* The last step is the terminal step. */
    *a_x = time_step[(istep >= nstep)?(nstep - 1):istep];
    *a_x1 = time_step[(istep + 1 >= nstep)?(nstep - 1):(istep + 1)];

    double a_xm1 = time_step[(istep > 0)?(istep - 1):0];
    *a_v = sqrt(a_xm1 * *(a_x));
    *a_v1 = sqrt(*a_x * *a_x1);

    VPM * vpm = vpm_find(fastpm->vpm_list, *a_x);
    fastpm->pm = &vpm->pm;
}

void 
fastpm_add_extension(FastPM * fastpm, 
    enum FastPMExtensionPoint where,
    void * function, void * userdata) 
{
    FastPMExtension * q = malloc(sizeof(FastPMExtension));
    q->userdata = userdata;
    q->function = function;
    q->next = fastpm->exts[where];
    fastpm->exts[where] = q;
}

