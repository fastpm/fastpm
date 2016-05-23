#include <string.h>
#include <alloca.h>
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
#include "solver-pm-internal.h"

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


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void
scale_acc(PMStore * po, double correction, double fudge);

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

    fastpm->model = malloc(sizeof(FastPMModel));

    if(!fastpm->USE_EXTERNAL_PSTORE) {
        fastpm->p = malloc(sizeof(PMStore));
        pm_store_init(fastpm->p);

        pm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3),
            PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC,
            fastpm->alloc_factor, comm);
    } else {
        fastpm->p = fastpm->USE_EXTERNAL_PSTORE;
    }

    fastpm->vpm_list = vpm_create(fastpm->vpminit,
                           &baseinit, &fastpm->p->iface, comm);

    fastpm_model_init(fastpm->model, fastpm, fastpm->USE_MODEL);

    int i = 0;
    for (i = 0; i < FASTPM_EXT_MAX; i ++) {
        fastpm->exts[i] = NULL;
    }

    PMInit pm_2lptinit = {
            .Nmesh = fastpm->nc,
            .BoxSize = fastpm->boxsize,
            .NprocY = 0, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = 0,
        };

    fastpm->pm_2lpt = malloc(sizeof(PM));
    pm_init(fastpm->pm_2lpt, &pm_2lptinit, &fastpm->p->iface, fastpm->comm);
}

void
fastpm_setup_ic(FastPM * fastpm, FastPMFloat * delta_k_ic)
{

    PM * pm_2lpt = fastpm->pm_2lpt;
    if(delta_k_ic) {
        double shift0;
        if(fastpm->USE_SHIFT) {
            shift0 = fastpm->boxsize / fastpm->nc * 0.5;
        } else {
            shift0 = 0;
        }
        double shift[3] = {shift0, shift0, shift0};

        int nc[3] = {fastpm->nc, fastpm->nc, fastpm->nc};

        pm_store_set_lagrangian_position(fastpm->p, pm_2lpt, shift, nc);

        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(pm_2lpt, delta_k_ic, fastpm->p, shift);
    }

    if(fastpm->USE_DX1_ONLY == 1) {
        memset(fastpm->p->dx2, 0, sizeof(fastpm->p->dx2[0]) * fastpm->p->np);
    }
    pm_store_summary(fastpm->p, pm_2lpt->Comm2D);
}

void
fastpm_evolve(FastPM * fastpm, double * time_step, int nstep) 
{
    CLOCK(decompose);
    CLOCK(force);
    CLOCK(kick);
    CLOCK(drift);
    CLOCK(afterforce);
    CLOCK(beforekick);
    CLOCK(beforedrift);
    CLOCK(correction);

    FastPMExtension * ext;

    fastpm_model_build(fastpm->model, fastpm->p, time_step[0], time_step[nstep - 1]);

    pm_2lpt_evolve(time_step[0], fastpm->p, fastpm->omega_m, fastpm->USE_DX1_ONLY);

    MPI_Barrier(fastpm->comm);

    CLOCK(warmup);

    if(fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(fastpm->p->v, 0, sizeof(fastpm->p->v[0]) * fastpm->p->np);
    }

    LEAVE(warmup);

    MPI_Barrier(fastpm->comm);

    double correction = 1.0;
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

        ENTER(force);
        fastpm_calculate_forces(fastpm, delta_k);
        LEAVE(force);

        ENTER(afterforce);
        for(ext = fastpm->exts[FASTPM_EXT_AFTER_FORCE];
            ext; ext = ext->next) {
                ((fastpm_ext_after_force) ext->function) 
                    (fastpm, delta_k, a_x, ext->userdata);
        }
        LEAVE(afterforce);


        pm_free(fastpm->pm, delta_k);

        /* correct for linear theory before kick and drift */
        ENTER(correction);
        fastpm_model_evolve(fastpm->model, a_x1);

        correction = fastpm_model_find_correction(fastpm->model, a_x, a_x1, a_v, a_v1);
        scale_acc(fastpm->p, correction, 1.0);
        /* the value of correction is leaked to the next step for reporting. */

        LEAVE(correction);

        ENTER(beforekick);
        FastPMKick kick;
        fastpm_kick_init(&kick, fastpm, fastpm->p, a_v1);
        FastPMDrift drift;
        fastpm_drift_init(&drift, fastpm, fastpm->p, a_x1);

        /* take snapshots if needed, before the kick */
        for(ext = fastpm->exts[FASTPM_EXT_BEFORE_KICK];
            ext; ext = ext->next) {
                ((fastpm_ext_before_kick) ext->function) 
                    (fastpm, &kick, ext->userdata);
        }
        LEAVE(beforekick);

        /* never go beyond 1.0 */
        if(a_x >= 1.0) break; 

        // Leap-frog "kick" -- velocities updated

        ENTER(kick);
        fastpm_kick_store(fastpm, fastpm->p, fastpm->p, a_v1);
        LEAVE(kick);

        ENTER(beforedrift);
        /* take snapshots if needed, before the drift */
        for(ext = fastpm->exts[FASTPM_EXT_BEFORE_DRIFT];
            ext; ext = ext->next) {
                ((fastpm_ext_before_drift) ext->function)
                    (fastpm, &drift, ext->userdata);
        }
        LEAVE(beforedrift);

        // Leap-frog "drift" -- positions updated

        ENTER(drift);
        fastpm_drift_store(fastpm, fastpm->p, fastpm->p, a_x1);
        LEAVE(drift);
    }

}

void
fastpm_destroy(FastPM * fastpm) 
{
    pm_destroy(fastpm->pm_2lpt);
    free(fastpm->pm_2lpt);
    fastpm_model_destroy(fastpm->model);
    if(!fastpm->USE_EXTERNAL_PSTORE)
        pm_store_destroy(fastpm->p);
    vpm_free(fastpm->vpm_list);
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
    /* skip scaling if there is no correction. */
    if(correction == 1.0) return;
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

