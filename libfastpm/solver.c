#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"

static void
fastpm_decompose(FastPMSolver * fastpm);


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

/* Useful stuff */
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);

void fastpm_solver_init(FastPMSolver * fastpm,
    FastPMConfig * config,
    MPI_Comm comm) {

    fastpm->config[0] = *config;

    fastpm->gravity[0] = (FastPMGravity) {
        .PainterType = config->PAINTER_TYPE,
        .PainterSupport = config->painter_support,
        .KernelType = config->KERNEL_TYPE,
        .DealiasingType = config->DEALIASING_TYPE,
        .ComputePotential = config->COMPUTE_POTENTIAL,
        .ComputeTidal = config->COMPUTE_TIDAL,
    };

    fastpm->cosmology[0] = (FastPMCosmology) {
        .OmegaM = config->omega_m,
        .OmegaLambda = 1.0 - config->omega_m,
    };

    fastpm->event_handlers = NULL;

    PMInit baseinit = {
            .Nmesh = config->nc,
            .BoxSize = config->boxsize,
            .NprocY = config->NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = config->UseFFTW,
        };

    fastpm->comm = comm;

    MPI_Comm_rank(comm, &fastpm->ThisTask);
    MPI_Comm_size(comm, &fastpm->NTask);

    fastpm->p = malloc(sizeof(FastPMStore));

    fastpm_store_init_evenly(fastpm->p, pow(1.0 * config->nc, 3),
          PACK_POS | PACK_VEL | PACK_ID
        | PACK_DX1 | PACK_DX2 | PACK_ACC
        | (config->SAVE_Q?PACK_Q:0)
        | (config->COMPUTE_POTENTIAL?PACK_POTENTIAL:0)
        | (config->COMPUTE_TIDAL?PACK_TIDAL:0),
        config->alloc_factor, comm);

    fastpm->vpm_list = vpm_create(config->vpminit,
                           &baseinit, comm);

    PMInit basepminit = {
            .Nmesh = config->nc,
            .BoxSize = config->boxsize,
            .NprocY = 0, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = 0,
        };

    fastpm->basepm = malloc(sizeof(PM));
    pm_init(fastpm->basepm, &basepminit, fastpm->comm);
}

void
fastpm_solver_setup_ic(FastPMSolver * fastpm, FastPMFloat * delta_k_ic)
{

    PM * basepm = fastpm->basepm;
    FastPMConfig * config = fastpm->config;
    FastPMStore * p = fastpm->p;

    if(delta_k_ic) {
        double shift0;
        if(config->USE_SHIFT) {
            shift0 = config->boxsize / config->nc * 0.5;
        } else {
            shift0 = 0;
        }
        double shift[3] = {shift0, shift0, shift0};

        int nc[3] = {config->nc, config->nc, config->nc};

        fastpm_store_set_lagrangian_position(p, basepm, shift, nc);

        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(basepm, delta_k_ic, p, shift);
    }

    if(config->USE_DX1_ONLY == 1) {
        memset(p->dx2, 0, sizeof(p->dx2[0]) * p->np);
    }
    fastpm_store_summary(p, fastpm->info.dx1, fastpm->info.dx2, fastpm->comm);
}
static void
fastpm_do_warmup(FastPMSolver * fastpm, double a0);
static void
fastpm_do_kick(FastPMSolver * fastpm, FastPMTransition * trans);
static void
fastpm_do_drift(FastPMSolver * fastpm, FastPMTransition * trans);
static void
fastpm_do_force(FastPMSolver * fastpm, FastPMTransition * trans);

static void
fastpm_do_interpolation(FastPMSolver * fastpm,
        FastPMDriftFactor * drift, FastPMKickFactor * kick, double a1, double a2);

void
fastpm_solver_evolve(FastPMSolver * fastpm, double * time_step, int nstep) 
{

    fastpm_do_warmup(fastpm, time_step[0]);

    FastPMStates * states = malloc(sizeof(FastPMStates));

    FastPMState template[] = {
    {0, 0, 1}, /* Kick */
    {0, 1, 1}, /* Drift */
    {0, 2, 1}, /* Drift */
    {2, 2, 1}, /* Force */
    {2, 2, 2}, /* Kick */
    {-1, -1, -1} /* End of table */
    };

    fastpm_tevo_generate_states(states, nstep-1, template, time_step);

    FastPMTransition transition[1];

    /* The last step is the 'terminal' step */
    int i;
    for(i = 1; states->table[i].force != -1; i ++) {
        fastpm_tevo_transition_init(transition, states, i - 1, i);

        FastPMTransitionEvent event[1];
        event->transition = transition;

        CLOCK(beforetransit);
        ENTER(beforetransit);
        fastpm_solver_emit_event(fastpm, FASTPM_EVENT_TRANSITION,
                FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) event);
        LEAVE(beforetransit);

        switch(transition->action) {
            case FASTPM_ACTION_KICK:
                fastpm_do_kick(fastpm, transition);
            break;
            case FASTPM_ACTION_DRIFT:
                fastpm_do_drift(fastpm, transition);
            break;
            case FASTPM_ACTION_FORCE:
                fastpm_do_force(fastpm, transition);
            break;
        }

        CLOCK(aftertransit);
        ENTER(aftertransit);
        fastpm_solver_emit_event(fastpm, FASTPM_EVENT_TRANSITION,
                FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event);
        LEAVE(aftertransit);

        if(i == 1) {
            /* Special treatment on the initial state because the
             * interpolation ranges are semi closed -- (, ] . we miss the initial step otherwise.
             * this needs to be after force calculation for potential to be valid. */
            double a0 = time_step[0];
            FastPMKickFactor kick;
            FastPMDriftFactor drift;
            fastpm_kick_init(&kick, fastpm, a0, a0, a0);
            fastpm_drift_init(&drift, fastpm, a0, a0, a0);
            fastpm_do_interpolation(fastpm, &drift, &kick, a0, a0);

        }
    }
    fastpm_tevo_destroy_states(states);
}

static void
fastpm_do_interpolation(FastPMSolver * fastpm,
        FastPMDriftFactor * drift, FastPMKickFactor * kick, double a1, double a2)
{
    CLOCK(interpolation);

    ENTER(interpolation);

    FastPMInterpolationEvent event[1];
    event->drift = drift;
    event->kick = kick;
    event->a1 = a1;
    event->a2 =a2;
    fastpm_solver_emit_event(fastpm, FASTPM_EVENT_INTERPOLATION, FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) event);

    LEAVE(interpolation);
}

static void
fastpm_do_warmup(FastPMSolver * fastpm, double a0)
{
    FastPMConfig * config = fastpm->config;
    CLOCK(warmup);
    ENTER(warmup);

    pm_2lpt_evolve(a0, fastpm->p, fastpm->cosmology, config->USE_DX1_ONLY);

    LEAVE(warmup);
}

PM *
fastpm_find_pm(FastPMSolver * fastpm, double a)
{
    VPM * vpm = vpm_find(fastpm->vpm_list, a);
    return &vpm->pm;
}

static void
fastpm_do_force(FastPMSolver * fastpm, FastPMTransition * trans)
{
    FastPMGravity * gravity = fastpm->gravity;

    CLOCK(decompose);
    CLOCK(force);
    CLOCK(afterforce);

    fastpm->pm = fastpm_find_pm(fastpm, trans->a.f);

    PM * pm = fastpm->pm;

    fastpm->info.Nmesh = fastpm->pm->init.Nmesh;

    FastPMFloat * delta_k = pm_alloc(pm);
    ENTER(decompose);
    fastpm_decompose(fastpm);
    LEAVE(decompose);

    ENTER(force);

    fastpm_gravity_calculate(gravity, pm, fastpm->p, delta_k);
    LEAVE(force);

    ENTER(afterforce);

    FastPMForceEvent event[1];
    event->delta_k = delta_k;
    event->a_f = trans->a.f;

    fastpm_solver_emit_event(fastpm, FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event);
    LEAVE(afterforce);

    pm_free(pm, delta_k);

}

static void
fastpm_do_kick(FastPMSolver * fastpm, FastPMTransition * trans)
{

    FastPMStore * p = fastpm->p;

    CLOCK(kick);

    FastPMKickFactor kick;
    fastpm_kick_init(&kick, fastpm, trans->a.i, trans->a.r, trans->a.f);
    if(trans->end->v == trans->end->x) {
        FastPMDriftFactor drift;

        FastPMTransition dual[1];
        fastpm_tevo_transition_find_dual(trans, dual);

        fastpm_drift_init(&drift, fastpm, dual->a.i, dual->a.r, dual->a.f);

        fastpm_do_interpolation(fastpm, &drift, &kick, trans->a.i, trans->a.f);
    }

    /* Do kick */
    ENTER(kick);
    fastpm_kick_store(&kick, p, p, trans->a.f);
    LEAVE(kick);
}

static void
fastpm_do_drift(FastPMSolver * fastpm, FastPMTransition * trans)
{
    FastPMStore * p = fastpm->p;

    CLOCK(drift);

    FastPMDriftFactor drift;
    fastpm_drift_init(&drift, fastpm, trans->a.i, trans->a.r, trans->a.f);

    if(trans->end->v == trans->end->x) {
        FastPMKickFactor kick;

        FastPMTransition dual[1];
        fastpm_tevo_transition_find_dual(trans, dual);
        fastpm_kick_init(&kick, fastpm, dual->a.i, dual->a.r, dual->a.f);

        fastpm_do_interpolation(fastpm, &drift, &kick, trans->a.i, trans->a.f);
    }

    /* Do drift */
    ENTER(drift);
    fastpm_drift_store(&drift, p, p, trans->a.f);
    LEAVE(drift);
}

void
fastpm_solver_destroy(FastPMSolver * fastpm) 
{
    pm_destroy(fastpm->basepm);
    free(fastpm->basepm);
    fastpm_store_destroy(fastpm->p);

    vpm_free(fastpm->vpm_list);

    /* FIXME: free VPM and stuff. */
    FastPMEventHandler * h, * h2;
    for(h = fastpm->event_handlers; h; h = h2) {
        h2 = h->next;
        free(h);
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

    int NTask;
    MPI_Comm_size(fastpm->comm, &NTask);

    /* apply periodic boundary and move particles to the correct rank */
    fastpm_store_wrap(fastpm->p, pm->BoxSize);
    fastpm_store_decompose(fastpm->p, to_rank, pm, fastpm->comm);
    size_t np_max;
    size_t np_min;


    double np_mean = 1.0 * fastpm_store_get_np_total(p, fastpm->comm) / NTask;

    MPI_Allreduce(&p->np, &np_max, 1, MPI_LONG, MPI_MAX, fastpm->comm);
    MPI_Allreduce(&p->np, &np_min, 1, MPI_LONG, MPI_MIN, fastpm->comm);

    fastpm->info.imbalance.min = np_min / np_mean;
    fastpm->info.imbalance.max = np_max / np_mean;
}

/* Interpolate position and velocity for snapshot at a=aout */
void
fastpm_set_snapshot(FastPMSolver * fastpm,
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                FastPMStore * po,
                double aout)
{
    FastPMStore * p = fastpm->p;
    FastPMCosmology * c = fastpm->cosmology;
    PM * pm = fastpm->basepm;
    int np = p->np;

    fastpm_kick_store(kick, p, po, aout);

    fastpm_drift_store(drift, p, po, aout);

    int i;
    /* potfactor converts fastpm Phi to dimensionless */
    double potfactor = 1.5 * c->OmegaM / (HubbleDistance * HubbleDistance);
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* convert the unit from a**2 dx/dt / H0 in Mpc/h to a dx/dt km/s */
            po->v[i][d] *= HubbleConstant / aout;
        }
        po->id[i] = p->id[i];
        /* convert the unit from comoving (Mpc/h) ** 2 to dimensionless potential. */
        if(po->potential)
            po->potential[i] = p->potential[i] / aout * potfactor;
        if(po->tidal) {
            for( d = 0; d < 3; d ++) {
                po->tidal[i][d] = p->tidal[i][d] / aout * potfactor;
            }
        }
    }

    po->np = np;
    po->a_x = po->a_v = aout;

    fastpm_store_wrap(po, pm->BoxSize);
    fastpm_store_decompose(po, to_rank, pm, fastpm->comm);
}

