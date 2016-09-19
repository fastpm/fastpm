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

    fastpm->p = malloc(sizeof(FastPMStore));
    fastpm_store_init(fastpm->p);

    fastpm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3),
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC | (fastpm->SAVE_Q?PACK_Q:0),
        fastpm->alloc_factor, comm);

    fastpm->vpm_list = vpm_create(fastpm->vpminit,
                           &baseinit, comm);

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

static void
fastpm_do_kick(FastPMSolver * fastpm, FastPMTETransition * thiskick);
static void
fastpm_do_drift(FastPMSolver * fastpm, FastPMTETransition * thisdrift);
static void
fastpm_do_force(FastPMSolver * fastpm, FastPMTETransition * thisforce);

void
fastpm_evolve(FastPMSolver * fastpm, double * time_step, int nstep) 
{
    MPI_Comm comm = fastpm->comm;

    pm_2lpt_evolve(time_step[0], fastpm->p, fastpm->omega_m, fastpm->USE_DX1_ONLY);

    MPI_Barrier(comm);

    CLOCK(warmup);

    if(fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(fastpm->p->v, 0, sizeof(fastpm->p->v[0]) * fastpm->p->np);
    }

    LEAVE(warmup);

    FastPMTEStates * states = malloc(sizeof(FastPMTEStates));

    FastPMTEEntry template[] = {
    {0, 0, 1}, /* Kick */
    {0, 1, 1}, /* Drift */
    {0, 2, 1}, /* Drift */
    {2, 2, 1}, /* Force */
    {2, 2, 2}, /* Kick */
    {-1, -1, -1} /* End of table */
    };

    fastpm_tevo_generate_states(states, nstep-1, template, time_step);

    FastPMTETransition thisdrift, thiskick, thisforce;
    FastPMTETransition lastdrift, lastkick;

    enum {
        FORCE = 1,
        KICK = 2,
        DRIFT = 3,
    } action;

    MPI_Barrier(comm);

    /* The last step is the 'terminal' step */
    int i = 1;
    while(states->table[i].a != -1) {
        if(states->table[i].a != states->table[i-1].a) {
            /* Force */
            action = FORCE;
            fastpm_tevo_transition_init(&thisforce, states,
                    states->table[i-1].a,
                    states->table[i].x,
                    states->table[i].a);

            fastpm->info.istep = i;
            fastpm->info.a_x = thisforce.a_i;
            fastpm->info.a_x1 = thisforce.a_f;
            fastpm->info.a_v =  thisforce.a_r;
        }
        if(states->table[i].v != states->table[i-1].v) {
            /* Kick */
            fastpm_tevo_transition_init(&thiskick, states,
                    states->table[i-1].v,
                    states->table[i-1].a,
                    states->table[i].v);

            action = KICK;
        }
        if(states->table[i].x != states->table[i-1].x) {
            /* Drift */
            fastpm_tevo_transition_init(&thisdrift, states,
                    states->table[i-1].x,
                    states->table[i-1].v,
                    states->table[i].x);
            action = DRIFT;
        }

        if(states->table[i].x == states->table[i].v && states->table[i].x != states->table[i].a) {
            /* Interpolation */
            if(action == KICK) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    thiskick.f, thiskick.i, thiskick.r,
                    lastdrift.i, lastdrift.f, lastdrift.r
                );
            }
            if(action == DRIFT) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    lastkick.i, lastkick.f, lastkick.r,
                    thisdrift.f, thisdrift.i, thisdrift.r
                );
            }
        }

        switch(action) {
            case FORCE:
                /* Calculate forces */
                fastpm_do_force(fastpm, &thisforce);
                break;

            case KICK:
                lastkick = thiskick;
                fastpm_do_kick(fastpm, &thiskick);
                break;

            case DRIFT:
                lastdrift = thisdrift;
                fastpm_do_drift(fastpm, &thisdrift);
                break;
        }
        i++;
    }

}
static void
fastpm_do_force(FastPMSolver * fastpm, FastPMTETransition * thisforce)
{

    FastPMExtension * ext;

    CLOCK(decompose);
    CLOCK(force);
    CLOCK(afterforce);

    VPM * vpm = vpm_find(fastpm->vpm_list, thisforce->a_r);
    fastpm->pm = &vpm->pm;
    PM * pm = fastpm->pm;

    fastpm->info.Nmesh = fastpm->pm->init.Nmesh;
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
            (fastpm, delta_k, thisforce->a_r, ext->userdata);
    }
    LEAVE(afterforce);

    pm_free(pm, delta_k);

}

static void
fastpm_do_kick(FastPMSolver * fastpm, FastPMTETransition * thiskick)
{
    FastPMExtension * ext;
    FastPMStore * p = fastpm->p;

    CLOCK(kick);
    CLOCK(beforekick);

    /* Used by callbacks */
    FastPMKick kick;
    fastpm_kick_init(&kick, fastpm, thiskick->a_i, thiskick->a_r, thiskick->a_f);

    ENTER(beforekick);
    /* take snapshots if needed, before the kick */
    for(ext = fastpm->exts[FASTPM_EXT_BEFORE_KICK];
            ext; ext = ext->next) {
        ((fastpm_ext_before_kick) ext->function) 
            (fastpm, &kick, ext->userdata);
    }
    LEAVE(beforekick);
    /* Do kick */
    ENTER(kick);
    fastpm_kick_store(fastpm, p, p, thiskick->a_f);
    LEAVE(kick);
}

static void
fastpm_do_drift(FastPMSolver * fastpm, FastPMTETransition * thisdrift)
{
    FastPMExtension * ext;
    FastPMStore * p = fastpm->p;

    CLOCK(drift);
    CLOCK(beforedrift);

    ENTER(beforedrift);

    /* Used by callbacks */
    FastPMDrift drift;
    fastpm_drift_init(&drift, fastpm, thisdrift->a_i, thisdrift->a_r, thisdrift->a_f);

    /* take snapshots if needed, before the drift */
    for(ext = fastpm->exts[FASTPM_EXT_BEFORE_DRIFT];
            ext; ext = ext->next) {
        ((fastpm_ext_before_drift) ext->function)
            (fastpm, &drift, ext->userdata);
    }
    LEAVE(beforedrift);
    /* Do drift */
    ENTER(drift);
    fastpm_drift_store(fastpm, p, p, thisdrift->a_f);
    LEAVE(drift);
}

void
fastpm_destroy(FastPMSolver * fastpm) 
{
    pm_destroy(fastpm->basepm);
    free(fastpm->basepm);
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
