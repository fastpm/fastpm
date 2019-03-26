#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/store.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"

static void
fastpm_decompose(FastPMSolver * fastpm);


#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

void fastpm_solver_init(FastPMSolver * fastpm,
    FastPMConfig * config,
    MPI_Comm comm) {

    fastpm->config[0] = *config;

    fastpm->cosmology[0] = (FastPMCosmology) {
        .OmegaM = config->omega_m,
        .OmegaLambda = 1.0 - config->omega_m,
    };

    if(config->pgdc)
    {
        fastpm->pgdc[0] = (FastPMPGDCorrection) {
	    .PainterType = config->PAINTER_TYPE,
	    .PainterSupport = config->painter_support,
	    .alpha0 = config->pgdc_alpha0,
	    .A = config->pgdc_A,
	    .B = config->pgdc_B,
	    .kl = config->pgdc_kl,
	    .ks = config->pgdc_ks,
        };
    }

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

    if(config->FORCE_TYPE == FASTPM_FORCE_COLA) {
        /* Cola requires DX1 and DX2 to be permantly stored. */
        config->ExtraAttributes |= COLUMN_DX1;
        config->ExtraAttributes |= COLUMN_DX2;
    }
    
    memset(fastpm->has_species, 0, FASTPM_SOLVER_NSPECIES);
    
    fastpm_store_init_evenly(fastpm->cdm,
          fastpm_species_get_name(FASTPM_SPECIES_CDM),
          pow(1.0 * config->nc, 3),
          COLUMN_POS | COLUMN_VEL | COLUMN_ID | COLUMN_MASK | COLUMN_ACC | fastpm->config->ExtraAttributes,
          fastpm->config->alloc_factor,
          fastpm->comm);

    fastpm_solver_add_species(fastpm, FASTPM_SPECIES_CDM, fastpm->cdm);   //add CDM [why make np_total a double?]
    
    fastpm->vpm_list = vpm_create(config->vpminit,
                           &baseinit, comm);

    PMInit basepminit = {
            .Nmesh = config->nc,
            .BoxSize = config->boxsize,
            .NprocY = config->NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 0, /* use untransposed to make sure we see all kz on a rank; this speeds up IC */
            .use_fftw = config->UseFFTW,
        };

    fastpm->basepm = malloc(sizeof(PM));
    pm_init(fastpm->basepm, &basepminit, fastpm->comm);
    
    double shift0;
    if(config->USE_SHIFT) {
        shift0 = config->boxsize / config->nc * 0.5;
    } else {
        shift0 = 0;
    }
    double shift[3] = {shift0, shift0, shift0};

    fastpm_store_fill(fastpm_solver_get_species(fastpm, FASTPM_SPECIES_CDM), fastpm->basepm, shift, NULL);
}

void
fastpm_solver_setup_lpt(FastPMSolver * fastpm,
        enum FastPMSpecies species,
        FastPMFloat * delta_k_ic,
        double a0)
{

    FastPMStore * p = fastpm_solver_get_species(fastpm, species);
    if(!p) fastpm_raise(-1, "Species requested (%d) does not exist", species);

    PM * basepm = fastpm->basepm;
    FastPMConfig * config = fastpm->config;

    double BoxSize = fastpm->config->boxsize;
    uint64_t NC = fastpm->config->nc;
    double OmegaM = fastpm->cosmology->OmegaM;
    double M0 = OmegaM * FASTPM_CRITICAL_DENSITY * (BoxSize / NC) * (BoxSize / NC) * (BoxSize / NC);

    p->meta.M0 = M0;

    int temp_dx1 = 0;
    int temp_dx2 = 0;
    if(p->dx1 == NULL) {
        p->dx1 = fastpm_memory_alloc(p->mem, "DX1", sizeof(p->dx1[0]) * p->np_upper, FASTPM_MEMORY_STACK);
        temp_dx1 = 1;
    }
    if(p->dx2 == NULL) {
        p->dx2 = fastpm_memory_alloc(p->mem, "DX2", sizeof(p->dx2[0]) * p->np_upper, FASTPM_MEMORY_STACK);
        temp_dx2 = 1;
    }

    FastPMLPTEvent event[1];
    event->pm = basepm;
    event->delta_k = delta_k_ic;
    event->p = p;

    fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_LPT,
                FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) event, fastpm);

    if(delta_k_ic) {
        double shift0;
        if(config->USE_SHIFT) {
            shift0 = config->boxsize / config->nc * 0.5;
        } else {
            shift0 = 0;
        }
        double shift[3] = {shift0, shift0, shift0};
        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(basepm, delta_k_ic, p, shift);
    }

    if(config->USE_DX1_ONLY == 1) {
        memset(p->dx2, 0, sizeof(p->dx2[0]) * p->np);
    }

    pm_2lpt_evolve(a0, p, fastpm->cosmology, config->USE_DX1_ONLY);

    fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_LPT,
                FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event, fastpm);

    if(temp_dx2) {
        fastpm_memory_free(p->mem, p->dx2);
        p->dx2 = NULL;
    }
    if(temp_dx1) {
        fastpm_memory_free(p->mem, p->dx1);
        p->dx1 = NULL;
    }
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

enum FastPMSpecies
fastpm_solver_iter_species(FastPMSolver * fastpm, int * iter)
{
    while(!fastpm->has_species[*iter]) {
        (*iter) ++;
        if(*iter >= sizeof(fastpm->has_species)) {
            return -1;
        }
    }
    enum FastPMSpecies sp = *iter;
    (*iter) ++;
    return sp;
}

FastPMStore *
fastpm_solver_get_species(FastPMSolver * fastpm, enum FastPMSpecies species)
{
    if(fastpm->has_species[species]) {
        return fastpm->species[species];
    } else {
        return NULL;
    }
}

void
fastpm_solver_add_species(FastPMSolver * fastpm, enum FastPMSpecies species, FastPMStore * store)
{   
    /*Adds a particle [store] of species type species to the solver. */

    fastpm->species[species] = store;
    fastpm->has_species[species] = 1;
}


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
        fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_TRANSITION,
                FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) event, fastpm);
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
        fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_TRANSITION,
                FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event, fastpm);
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
    free(states);
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
    fastpm_emit_event(fastpm->event_handlers,
            FASTPM_EVENT_INTERPOLATION, FASTPM_EVENT_STAGE_BEFORE,
            (FastPMEvent*) event, fastpm);

    LEAVE(interpolation);
}

static void
fastpm_do_warmup(FastPMSolver * fastpm, double a0)
{
    CLOCK(warmup);
    ENTER(warmup);

    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;
        /* set acc to zero or we see valgrind errors */
        memset(p->acc, 0, sizeof(p->acc[0]) * p->np);
    }
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
    CLOCK(decompose);
    CLOCK(force);
    CLOCK(pgdc);
    CLOCK(event);

    fastpm->pm = fastpm_find_pm(fastpm, trans->a.f);

    FastPMPainter painter[1];
    PM * pm = fastpm->pm;

    FastPMFloat * delta_k = pm_alloc(pm);

    FastPMForceEvent event[1];

    /* FIXME modify gravity.c to compute delta_k from all species. */
    FastPMStore * p = fastpm_solver_get_species(fastpm, FASTPM_SPECIES_CDM);

    fastpm_painter_init(painter, pm, fastpm->config->PAINTER_TYPE, fastpm->config->painter_support);

    int64_t N = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_LONG, MPI_SUM, fastpm->comm);

    event->delta_k = delta_k;
    event->a_f = trans->a.f;
    event->pm = pm;
    event->N = N; /* FIXME: pass in the shot noise instead? */
    event->painter = painter;
    event->kernel = fastpm->config->KERNEL_TYPE;

    /* find the time stamp of the next force calculation. This will
     * be useful for interpolating potentials of the structured mesh */
    FastPMTransition next[1];

    if(!fastpm_tevo_transition_find_next(trans, next)) {
        event->a_n = -1;
    } else {
        if(next->a.i != trans->a.f) {
            fastpm_raise(-1, "Failed to find next Force calculation\n");
        }
        event->a_n = next->a.f;
    }

    ENTER(decompose);
    fastpm_decompose(fastpm);
    LEAVE(decompose);

    fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) event, fastpm);

    ENTER(force);
    fastpm_solver_compute_force(fastpm, painter, fastpm->config->DEALIASING_TYPE, fastpm->config->KERNEL_TYPE, delta_k);
    LEAVE(force);

    if(p->pgdc) {
        ENTER(pgdc);
        /* delta_k is input; unchanged */
        FastPMPGDCorrection * pgdc = fastpm->pgdc;
        fastpm_pgdc_calculate(pgdc, pm, p, delta_k, trans->a.f, 1.0);
        LEAVE(pgdc);
    }
    ENTER(event);

    /* Prepare for the event; apply CIC correction for the painting, such that
     * the event has a compensated power spectrum. See e.g. MP-Gadget's gravpm.c */

    /* FIXME: the compensation transfer shall be a method of the painter. */
    fastpm_apply_decic_transfer(pm, delta_k, delta_k);

    fastpm_emit_event(fastpm->event_handlers, FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event, fastpm);
    LEAVE(event);

    pm_free(pm, delta_k);

}

static void
fastpm_do_kick(FastPMSolver * fastpm, FastPMTransition * trans)
{

    CLOCK(kick);

    FastPMKickFactor kick;
    fastpm_kick_init(&kick, fastpm, trans->a.i, trans->a.r, trans->a.f);
    if(trans->end->v == trans->end->x) {
        FastPMDriftFactor drift;

        FastPMTransition dual[1];
        if(!fastpm_tevo_transition_find_dual(trans, dual)) {
            fastpm_raise(-1, "Dual transition not found. The state table is likely wrong. Look at states->table.\n");
        }

        fastpm_drift_init(&drift, fastpm, dual->a.i, dual->a.r, dual->a.f);

        fastpm_do_interpolation(fastpm, &drift, &kick, trans->a.i, trans->a.f);
    }

    /* Do kick */
    ENTER(kick);
    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;
        
        if(kick.ai != p->meta.a_v) {
            fastpm_raise(-1, "kick is inconsitant with state.\n");
        }
        if(kick.ac != p->meta.a_x) {
            fastpm_raise(-1, "kick is inconsitant with state.\n");
        }
        fastpm_kick_store(&kick, p, p, trans->a.f);
    }
    LEAVE(kick);
}

static void
fastpm_do_drift(FastPMSolver * fastpm, FastPMTransition * trans)
{
    CLOCK(drift);

    FastPMDriftFactor drift;
    fastpm_drift_init(&drift, fastpm, trans->a.i, trans->a.r, trans->a.f);

    if(trans->end->v == trans->end->x) {
        FastPMKickFactor kick;

        FastPMTransition dual[1];
        if(!fastpm_tevo_transition_find_dual(trans, dual)) {
            fastpm_raise(-1, "Dual transition not found. The state table is likely wrong. Look at states->table.\n");
        }
        fastpm_kick_init(&kick, fastpm, dual->a.i, dual->a.r, dual->a.f);

        fastpm_do_interpolation(fastpm, &drift, &kick, trans->a.i, trans->a.f);
    }

    /* Do drift */
    ENTER(drift);
    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;

        if(drift.ai != p->meta.a_x) {
            fastpm_raise(-1, "drift is inconsitant with state.\n");
        }
        if(drift.ac != p->meta.a_v) {
            fastpm_raise(-1, "drift is inconsitant with state.\n");
        }
        fastpm_drift_store(&drift, p, p, trans->a.f);
        LEAVE(drift);
    }
}

void
fastpm_solver_destroy(FastPMSolver * fastpm) 
{
    pm_destroy(fastpm->basepm);
    free(fastpm->basepm);
    fastpm_store_destroy(fastpm->cdm);
    vpm_free(fastpm->vpm_list);

    fastpm_destroy_event_handlers(&fastpm->event_handlers);
}

static void
fastpm_decompose(FastPMSolver * fastpm) {
    PM * pm = fastpm->pm;

    int NTask;
    MPI_Comm_size(fastpm->comm, &NTask);

    int si;
    for(si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
        FastPMStore * p = fastpm_solver_get_species(fastpm, si);
        if(!p) continue;

        /* apply periodic boundary and move particles to the correct rank */
        fastpm_store_wrap(p, pm->BoxSize);

        if(0 != fastpm_store_decompose(p,
                (fastpm_store_target_func) FastPMTargetPM, pm,
                fastpm->comm))
        {
            fastpm_raise(-1, "Out of particle storage space\n");
        }
    }
}

/* Interpolate position and velocity for snapshot at a=aout,
 * this alters fastpm->species[], thus need to call unset to revert it.
 * 
 *
 * fastpm_set_snapshot(fastpm, .......)
 *
 *  Avoid using fastpm->species[] between
 *
 * fastpm_unset_snapshot(fastpm, .......)
 * */

void
fastpm_set_snapshot(FastPMSolver * fastpm,
                FastPMSolver * snapshot,
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                double aout) {

    int si;
    for (si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
        FastPMStore * p, *po;

        p  = fastpm_solver_get_species(fastpm, si);
        po = fastpm_solver_get_species(snapshot, si);

        if(!p || !po) { continue; }

        fastpm_set_species_snapshot(fastpm, p, drift, kick, po, aout);
    }
}

void
fastpm_unset_snapshot(FastPMSolver * fastpm,
                FastPMSolver * snapshot,
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                double aout) {

    int si;
    for (si = 0; si < FASTPM_SOLVER_NSPECIES; si ++) {
        FastPMStore * p, *po;

        p  = fastpm_solver_get_species(fastpm, si);
        po = fastpm_solver_get_species(snapshot, si);

        if(!p || !po) { continue; }

        fastpm_unset_species_snapshot(fastpm, p, drift, kick, po, aout);
    }
}


void
fastpm_set_species_snapshot(FastPMSolver * fastpm,
                FastPMStore * p,
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                FastPMStore * po,
                double aout)
{
    FastPMCosmology * c = fastpm->cosmology;
    PM * pm = fastpm->basepm;
    int np = p->np;

    memcpy(po, p, sizeof(FastPMStore));

    /* steal columns, but velocity, since we need to use old velocity
     * during the drift */
    fastpm_store_steal(p, po, p->attributes);

    /* Fake the attributes */
    po->attributes = p->attributes;

    /* update position; before kick to use the old velocity */
    fastpm_drift_store(drift, p, po, aout);

    /* update velocity */
    fastpm_kick_store(kick, p, po, aout);

    ptrdiff_t i;

    /* convert units */

    /* potfactor converts fastpm Phi to dimensionless */
    double potfactor = 1.5 * c->OmegaM / (HubbleDistance * HubbleDistance);

#pragma omp parallel for
    for(i=0; i<np; i++) {

        int d;
        for(d = 0; d < 3; d ++) {
            /* convert the unit from a**2 dx/dt / H0 in Mpc/h to a dx/dt km/s */
            po->v[i][d] *= HubbleConstant / aout;
        }

        /* convert the unit from comoving (Mpc/h) ** 2 to dimensionless potential. */
        if(po->potential)
            po->potential[i] *= potfactor / aout;
        if(po->tidal) {
            for( d = 0; d < 3; d ++) {
                po->tidal[i][d] *= potfactor / aout;
            }
        }
    }
    fastpm_store_wrap(po, pm->BoxSize);
}

/* revert the effect of a snapshot on fastpm->species, and destroy po */
void
fastpm_unset_species_snapshot(FastPMSolver * fastpm,
                FastPMStore * p,
                FastPMDriftFactor * drift,
                FastPMKickFactor * kick,
                FastPMStore * po,
                double aout)
{
    FastPMCosmology * c = fastpm->cosmology;
    PM * pm = fastpm->basepm;
    int np = po->np;

    ptrdiff_t i;

    /* convert units */

    /* potfactor converts fastpm Phi to dimensionless */
    double potfactor = 1.5 * c->OmegaM / (HubbleDistance * HubbleDistance);

#pragma omp parallel for
    for(i=0; i<np; i++) {

        int d;
        for(d = 0; d < 3; d ++) {
            /* convert the unit from a**2 dx/dt / H0 in Mpc/h to a dx/dt km/s */
            po->v[i][d] /= HubbleConstant / aout;
        }

        /* convert the unit from comoving (Mpc/h) ** 2 to dimensionless potential. */
        if(po->potential)
            po->potential[i] /= potfactor / aout;
        if(po->tidal) {
            for( d = 0; d < 3; d ++) {
                po->tidal[i][d] /= potfactor / aout;
            }
        }
    }

    /* revert velocity */
    fastpm_kick_store(kick, po, po, p->meta.a_v);

    /* revert position */
    fastpm_drift_store(drift, po, po, p->meta.a_x);

    /* steal back columns */
    fastpm_store_steal(po, p, p->attributes);

    fastpm_store_wrap(p, pm->BoxSize);

    /* Stop faking the attributes */
    po->attributes = 0;

}

double
fastpm_solver_growth_factor(FastPMSolver * fastpm, double a)
{
    return GrowthFactor(a, fastpm->cosmology);
}

double
fastpm_solver_growth_rate(FastPMSolver * fastpm, double a)
{
    return DLogGrowthFactor(a, fastpm->cosmology);
}
