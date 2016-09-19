#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/timemachine.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"
#include "solver-pm-internal.h"

static void
fastpm_decompose(FastPMSolver * fastpm);


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

/* Useful stuff */
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);

void fastpm_init(FastPMSolver * fastpm, 
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

//    fastpm->model = malloc(sizeof(FastPMModel));

    fastpm->p = malloc(sizeof(FastPMStore));
    fastpm_store_init(fastpm->p);

    fastpm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3),
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC | (fastpm->SAVE_Q?PACK_Q:0),
        fastpm->alloc_factor, comm);

    fastpm->vpm_list = vpm_create(fastpm->vpminit,
                           &baseinit, comm);

//    fastpm_model_init(fastpm->model, fastpm, fastpm->USE_MODEL);

    int i = 0;
    for (i = 0; i < FASTPM_EXT_MAX; i ++) {
        fastpm->exts[i] = NULL;
    }

    PMInit basepminit = {
            .Nmesh = fastpm->nc,
            .BoxSize = fastpm->boxsize,
            .NprocY = 0, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = 0,
        };

    fastpm->basepm = malloc(sizeof(PM));
    pm_init(fastpm->basepm, &basepminit, fastpm->comm);
}

void
fastpm_setup_ic(FastPMSolver * fastpm, FastPMFloat * delta_k_ic)
{

    PM * basepm = fastpm->basepm;
    FastPMStore * p = fastpm->p;

    if(delta_k_ic) {
        double shift0;
        if(fastpm->USE_SHIFT) {
            shift0 = fastpm->boxsize / fastpm->nc * 0.5;
        } else {
            shift0 = 0;
        }
        double shift[3] = {shift0, shift0, shift0};

        int nc[3] = {fastpm->nc, fastpm->nc, fastpm->nc};

        fastpm_store_set_lagrangian_position(p, basepm, shift, nc);

        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(basepm, delta_k_ic, p, shift);
    }

    if(fastpm->USE_DX1_ONLY == 1) {
        memset(p->dx2, 0, sizeof(p->dx2[0]) * p->np);
    }
    fastpm_store_summary(p, fastpm->info.dx1, fastpm->info.dx2, fastpm->comm);
}

void
fastpm_evolve(FastPMSolver * fastpm, double * time_step, int nstep) 
{
    CLOCK(decompose);
    CLOCK(force);
    CLOCK(kick);
    CLOCK(drift);
    CLOCK(afterforce);
    CLOCK(beforekick);
    CLOCK(beforedrift);

    FastPMExtension * ext;
    FastPMStore * p = fastpm->p;
    MPI_Comm comm = fastpm->comm;

//    fastpm_model_build(fastpm->model, time_step[0], time_step[nstep - 1]);

    pm_2lpt_evolve(time_step[0], p, fastpm->omega_m, fastpm->USE_DX1_ONLY);

    MPI_Barrier(comm);

    CLOCK(warmup);

    if(fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(p->v, 0, sizeof(p->v[0]) * p->np);
    }

    LEAVE(warmup);

    FastPMTEStates *states = malloc(sizeof(FastPMTEStates));

    FastPMTEEntry template[] = {
    {0, 0, 1}, /* Kick */
    {0, 1, 1}, /* Drift */
    {0, 2, 1}, /* Drift */
    {2, 2, 1}, /* Force */
    {2, 2, 2}, /* Kick */
    {-1, -1, -1} /* End of table */
    };

    fastpm_tevo_generate_states(states, nstep-1, template, time_step);

    fastpm_tevo_print_states(states);
//    FastPMTEStep lastkick = {-1, -1, -1}, lastdrift = {-1, -1, -1};
    FastPMTEStep thisdrift, thiskick;

    enum {
        FORCE = 1,
        KICK = 2,
        DRIFT = 3,
    } action;

    MPI_Barrier(comm);

    /* The last step is the 'terminal' step */
    int i = 1;
    while(states->table[i].a != -1) {
        fastpm->info.istep = i;
        fastpm->info.a_x = fastpm_tevo_i2t(states, states->table[i-1].x);
        fastpm->info.a_x1 = fastpm_tevo_i2t(states, states->table[i].x);
        fastpm->info.a_v =  fastpm_tevo_i2t(states, states->table[i-1].v);
        fastpm->info.a_v1 = fastpm_tevo_i2t(states, states->table[i].v);
    
        if(states->table[i].a != states->table[i-1].a) {
            /* Force */
            action = FORCE;
        }
        if(states->table[i].v != states->table[i-1].v) {
            /* Kick */
            thiskick.i = states->table[i-1].v;
            thiskick.f = states->table[i].v;
            thiskick.r = states->table[i-1].a;
            action = KICK;
        }
        if(states->table[i].x != states->table[i-1].x) {
            /* Drift */
            thisdrift.i = states->table[i-1].x;
            thisdrift.f = states->table[i].x;
            thisdrift.r = states->table[i-1].v;
            action = DRIFT;
        }

        switch(action) {
            case FORCE:
                /* Calculate forces */

                {
                double a_x = fastpm_tevo_i2t(states, states->table[i].x);

                VPM * vpm = vpm_find(fastpm->vpm_list, a_x);
                fastpm->pm = &vpm->pm;
                PM * pm = fastpm->pm;

                fastpm->info.Nmesh = fastpm->pm->init.Nmesh;
                fastpm_painter_init(fastpm->painter, pm, fastpm->PAINTER_TYPE, fastpm->painter_support);

                fastpm_painter_init(fastpm->painter, pm, fastpm->PAINTER_TYPE, fastpm->painter_support);

                FastPMFloat * delta_k = pm_alloc(pm);
                ENTER(decompose);
                fastpm_decompose(fastpm);
                LEAVE(decompose);

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

                pm_free(pm, delta_k);

                break;
                }
            case KICK:
                //lastkick = thiskick;

                ENTER(beforekick);

                /* Used by callbacks */
                FastPMKick kick;
                fastpm_kick_init(&kick, fastpm,
                        fastpm_tevo_i2t(states, thiskick.i),
                        fastpm_tevo_i2t(states, thiskick.r),
                        fastpm_tevo_i2t(states, thiskick.f));

                /* take snapshots if needed, before the kick */
                for(ext = fastpm->exts[FASTPM_EXT_BEFORE_KICK];
                    ext; ext = ext->next) {
                        ((fastpm_ext_before_kick) ext->function) 
                            (fastpm, &kick, ext->userdata);
                }
                LEAVE(beforekick);
                /* Do kick */
                ENTER(kick);
                fastpm_kick_store(fastpm, p, p, fastpm_tevo_i2t(states, thiskick.f));
                LEAVE(kick);
                break;

            case DRIFT:
                //lastdrift = thisdrift;

                ENTER(beforedrift);

                /* Used by callbacks */
                FastPMDrift drift;
                fastpm_drift_init(&drift, fastpm,
                        fastpm_tevo_i2t(states, thisdrift.i),
                        fastpm_tevo_i2t(states, thisdrift.r),
                        fastpm_tevo_i2t(states, thisdrift.f));

                /* take snapshots if needed, before the drift */
                for(ext = fastpm->exts[FASTPM_EXT_BEFORE_DRIFT];
                    ext; ext = ext->next) {
                        ((fastpm_ext_before_drift) ext->function)
                            (fastpm, &drift, ext->userdata);
                }
                LEAVE(beforedrift);
                /* Do drift */
                ENTER(drift);
                fastpm_drift_store(fastpm, p, p, fastpm_tevo_i2t(states, thisdrift.f));
                LEAVE(drift);

                break;
        }
        i++;
    }

}

void
fastpm_destroy(FastPMSolver * fastpm) 
{
    pm_destroy(fastpm->basepm);
    free(fastpm->basepm);
//    fastpm_model_destroy(fastpm->model);
    fastpm_store_destroy(fastpm->p);

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
    FastPMStore * p = (FastPMStore *) pdata;
    PM * pm = (PM*) data;
    double pos[3];
    p->get_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
}

static void
fastpm_decompose(FastPMSolver * fastpm) {
    PM * pm = fastpm->pm;
    FastPMStore * p = fastpm->p;
    /* apply periodic boundary and move particles to the correct rank */
    fastpm_store_wrap(fastpm->p, pm->BoxSize);
    fastpm_store_decompose(fastpm->p, to_rank, pm, fastpm->comm);
    size_t np_max;
    size_t np_min;

    /* FIXME move NTask to somewhere else. */
    double np_mean = pow(fastpm->nc, 3) / pm->NTask;
    MPI_Allreduce(&p->np, &np_max, 1, MPI_LONG, MPI_MAX, fastpm->comm);
    MPI_Allreduce(&p->np, &np_min, 1, MPI_LONG, MPI_MIN, fastpm->comm);

    fastpm->info.imbalance.min = np_min / np_mean;
    fastpm->info.imbalance.max = np_max / np_mean;
}

void 
fastpm_interp(FastPMSolver * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata) 
{
    /* interpolate and write snapshots, assuming p 
     * is at time a_x and a_v. */
    FastPMStore * p = fastpm->p;
    double a_x = p->a_x;
    double a_v = p->a_v;
    int iout;
    for(iout = 0; iout < nout; iout ++) {
        if(
        ! /* after a kick */
        (a_x < aout[iout] && aout[iout] < a_v)
        &&
        ! /* after a drift */
        (a_x >= aout[iout] && aout[iout] >= a_v)
        ) continue;

        FastPMStore * snapshot = alloca(sizeof(FastPMStore));
        fastpm_store_init(snapshot);
        fastpm_store_alloc(snapshot, p->np_upper, PACK_ID | PACK_POS | PACK_VEL);

        fastpm_info("Taking a snapshot...\n");

        fastpm_set_snapshot(fastpm, p, snapshot, aout[iout]);

        action(fastpm, snapshot, aout[iout], userdata);

        fastpm_store_destroy(snapshot);

    }
}

void 
fastpm_add_extension(FastPMSolver * fastpm, 
    enum FastPMExtensionPoint where,
    void * function, void * userdata) 
{
    FastPMExtension * q = malloc(sizeof(FastPMExtension));
    q->userdata = userdata;
    q->function = function;
    q->next = fastpm->exts[where];
    fastpm->exts[where] = q;
}
