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
fastpm_solver_setup_ic(FastPMSolver * fastpm, FastPMFloat * delta_k_ic)
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

        CLOCK(beforetransit);
        ENTER(beforetransit);
        FastPMExtension * ext;

        for(ext = fastpm->exts[FASTPM_EXT_BEFORE_TRANSITION];
                ext; ext = ext->next) {
            ((fastpm_ext_transition) ext->function)
                (fastpm, transition, ext->userdata);
        }
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
    }
    fastpm_tevo_destroy_states(states);
}

static void
fastpm_do_interpolation(FastPMSolver * fastpm,
        FastPMDriftFactor * drift, FastPMKickFactor * kick, double a1, double a2)
{
    CLOCK(interpolation);

    ENTER(interpolation);

    FastPMExtension * ext;
    for(ext = fastpm->exts[FASTPM_EXT_INTERPOLATE];
            ext; ext = ext->next) {
        ((fastpm_ext_interpolate) ext->function)
            (fastpm, drift, kick, a1, a2, ext->userdata);
    }
    LEAVE(interpolation);
}

static void
fastpm_do_warmup(FastPMSolver * fastpm, double a0)
{
    CLOCK(warmup);
    ENTER(warmup);

    pm_2lpt_evolve(a0, fastpm->p, fastpm->omega_m, fastpm->USE_DX1_ONLY);

    LEAVE(warmup);

    FastPMKickFactor kick;
    FastPMDriftFactor drift;
    fastpm_kick_init(&kick, fastpm, a0, a0, a0);
    fastpm_drift_init(&drift, fastpm, a0, a0, a0);
    fastpm_do_interpolation(fastpm, &drift, &kick, a0, a0);
}


static void
fastpm_do_force(FastPMSolver * fastpm, FastPMTransition * trans)
{

    FastPMGravity gravity = {
        .PainterType = fastpm->PAINTER_TYPE,
        .PainterSupport = fastpm->painter_support,
        .KernelType = fastpm->KERNEL_TYPE,
        .DealiasingType = fastpm->DEALIASING_TYPE,
    };

    FastPMExtension * ext;

    CLOCK(decompose);
    CLOCK(force);
    CLOCK(afterforce);

    VPM * vpm = vpm_find(fastpm->vpm_list, trans->a.f);
    fastpm->pm = &vpm->pm;
    PM * pm = fastpm->pm;

    fastpm->info.Nmesh = fastpm->pm->init.Nmesh;

    FastPMFloat * delta_k = pm_alloc(pm);
    ENTER(decompose);
    fastpm_decompose(fastpm);
    LEAVE(decompose);

    ENTER(force);

    fastpm_gravity_calculate(&gravity, pm, fastpm->p, delta_k);
    LEAVE(force);

    ENTER(afterforce);
    for(ext = fastpm->exts[FASTPM_EXT_AFTER_FORCE];
            ext; ext = ext->next) {
        ((fastpm_ext_after_force) ext->function)
            (fastpm, delta_k, trans->a.f, ext->userdata);
    }
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
    FastPMExtension * ext, * e2;
    int i;
    for(i = 0; i < FASTPM_EXT_MAX; i++) {
        for(ext = fastpm->exts[i]; ext; ext = e2) {
            e2 = ext->next;
            free(ext);
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
fastpm_solver_add_extension(FastPMSolver * fastpm,
    enum FastPMExtensionPoint where,
    void * function, void * userdata)
{
    FastPMExtension * q = malloc(sizeof(FastPMExtension));
    q->userdata = userdata;
    q->function = function;
    q->next = fastpm->exts[where];
    fastpm->exts[where] = q;
}

/* Interpolate position and velocity for snapshot at a=aout */
void
fastpm_set_snapshot(
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                FastPMStore * p, FastPMStore * po,
                double aout)
{
    int np= p->np;

    double H0 = 100.0f; // H0= 100 km/s/(h^-1 Mpc)

    fastpm_kick_store(kick, p, po, aout);

    fastpm_drift_store(drift, p, po, aout);

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* convert the unit from a**2 H_0 dx/dt in Mpc/h to a dx/dt km/s */
            po->v[i][d] *= H0 / aout;
        }
        po->id[i] = p->id[i];
    }

    po->np = np;
    po->a_x = po->a_v = aout;
}

