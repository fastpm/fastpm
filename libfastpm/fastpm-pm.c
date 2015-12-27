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
              double af);

void 
fastpm_drift(FastPM * fastpm,
               PMStore * pi, PMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPM * fastpm,
                PMStore * p, PMStore * po,
                double aout);

double fastpm_growth_factor(FastPM * fastpm, double a);


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void
fix_linear_growth(PMStore * p, double correction, double fudge);

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

    fastpm->comm = comm;
    MPI_Comm_rank(comm, &fastpm->ThisTask);
    MPI_Comm_size(comm, &fastpm->NTask);

    fastpm->p= malloc(sizeof(PMStore));
    fastpm->pm_2lpt = malloc(sizeof(PM));

    pm_store_init(fastpm->p);

    pm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC, 
        fastpm->alloc_factor, comm);

    fastpm->vpm_list = vpm_create(fastpm->vpminit,
                           &baseinit, &fastpm->p->iface, comm);

    pm_init_simple(fastpm->pm_2lpt, fastpm->p, fastpm->nc, fastpm->boxsize, comm);

    int i = 0;
    for (i = 0; i < FASTPM_EXT_MAX; i ++) {
        fastpm->exts[i] = NULL;
    }

}

void 
fastpm_solve_2lpt(FastPM * fastpm, FastPMFloat * delta_k_ic)
{
    double shift[3] = {
        fastpm->boxsize / fastpm->nc * 0.5,
        fastpm->boxsize / fastpm->nc * 0.5,
        fastpm->boxsize / fastpm->nc * 0.5,
        };

    pm_store_set_lagrangian_position(fastpm->p, fastpm->pm_2lpt, shift);

    /* read out values at locations with an inverted shift */
    pm_2lpt_solve(fastpm->pm_2lpt, delta_k_ic, fastpm->p, shift);
}

void
fastpm_evolve(FastPM * fastpm, double * time_step, int nstep) 
{
    FastPMExtension * ext;

    MPI_Barrier(fastpm->comm);

    CLOCK(warmup);

    pm_store_summary(fastpm->p, fastpm->comm);
    pm_2lpt_evolve(time_step[0], fastpm->p, fastpm->omega_m);
    if(fastpm->USE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(fastpm->p->v, 0, sizeof(fastpm->p->v[0]) * fastpm->p->np);
    }

    LEAVE(warmup);

    MPI_Barrier(fastpm->comm);

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

        CLOCK(decompose);
        fastpm_decompose(fastpm);

        LEAVE(decompose);

        /* Calculate PM forces. */
        FastPMFloat * delta_k = pm_alloc(fastpm->pm);

        /* watch out: boost the density since mesh is finer than grid */
        double density_factor = fastpm->pm->Norm / pow(1.0 * fastpm->nc, 3);

        CLOCK(force);
        pm_calculate_forces(fastpm->p, fastpm->pm, delta_k, density_factor);
        LEAVE(force);

        CLOCK(afterforce);
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_FORCE];
            ext; ext = ext->next) {
                ((fastpm_ext_after_force) ext->function) 
                    (fastpm, delta_k, a_x, ext->userdata);
        }
        LEAVE(afterforce);

        CLOCK(correction);
        double Plin = pm_calculate_linear_power(fastpm->pm, delta_k, fastpm->K_LINEAR);

        Plin /= pow(fastpm_growth_factor(fastpm, a_x), 2.0);
        if(istep == 0) {
            Plin0 = Plin;
        }

        double correction = sqrt(Plin0 / Plin);

        if(!fastpm->USE_LINEAR_THEORY) correction = 1.0;
        fastpm_info("<P(k<%g)> = %g Linear Theory = %g, correction=%g\n", 
                          fastpm->K_LINEAR, Plin, Plin0, correction); 
        fix_linear_growth(fastpm->p, correction, 2.0);
        LEAVE(correction);

        pm_free(fastpm->pm, delta_k);

        CLOCK(afterdrift);
        /* take snapshots if needed, before the kick */
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_DRIFT];
            ext; ext = ext->next) {
                ((fastpm_ext_after_drift) ext->function) 
                    (fastpm, ext->userdata);
        }
        LEAVE(afterdrift);

        /* never go beyond 1.0 */
        if(a_x >= 1.0) break; 
        
        // Leap-frog "kick" -- velocities updated

        CLOCK(kick);
        fastpm_kick(fastpm, fastpm->p, fastpm->p, a_v1);
        LEAVE(kick);

        CLOCK(afterkick);
        /* take snapshots if needed, before the drift */
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_KICK];
            ext; ext = ext->next) {
                ((fastpm_ext_after_kick) ext->function) 
                    (fastpm, ext->userdata);
        }
        LEAVE(afterkick);
        
        // Leap-frog "drift" -- positions updated

        CLOCK(drift);
        fastpm_drift(fastpm, fastpm->p, fastpm->p, a_x1);
        LEAVE(drift);

        /* no need to check for snapshots here, it will be checked next loop.  */
    }

}
void 
fastpm_destroy(FastPM * fastpm) 
{
    pm_store_destroy(fastpm->p);
    vpm_free(fastpm->vpm_list);
    pm_destroy(fastpm->pm_2lpt);
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
fix_linear_growth(PMStore * p, double correction, double fudge)
{
    ptrdiff_t i;
    int d;
    correction = pow(correction, fudge);
    for(d = 0; d < 3; d ++) {
        for(i = 0; i < p->np; i ++) {
            p->acc[i][d] *= correction;
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

