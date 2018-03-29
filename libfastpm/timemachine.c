#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>

static int
fastpm_tevo_block_len(FastPMState * template)
{
    /* Return number of elements in a template */
    int i = 0;
    while(template[i].force != -1) { i++; }
    return i++;
}

FastPMStates *
fastpm_tevo_generate_states(FastPMStates *states, int cycles, FastPMState *template, double *ts)
{
    /* Generate state table */
    int i, j, len = fastpm_tevo_block_len(template), N = len * cycles;

    double *timesteps = malloc((cycles + 1) * sizeof(double));
    FastPMState *table = malloc((N + 3) * sizeof(FastPMState));

    table[0].x = 0; // Initial conditions
    table[0].v = 0;
    table[0].force = -2;

    table[1].x = 0; // Force calculation
    table[1].v = 0;
    table[1].force = 0;

    for(i = 0; i < cycles; i++) { // Templates
        for(j = 0; j < len; j++) {
            table[j + i * len + 2].force = table[i * len + 1].force + template[j].force;
            table[j + i * len + 2].x = table[i * len + 1].x + template[j].x;
            table[j + i * len + 2].v = table[i * len + 1].v + template[j].v;
        }
    }

    table[N+2].force = -1; // End of table
    table[N+2].x = -1;
    table[N+2].v = -1;

    states->table = table; // Pack in struct
    states->cycle_len = template[len-1].force;
    states->cycles = cycles;

    for(i = 0; i <= cycles; i++) { timesteps[i] = ts[i]; }
    states->timesteps = timesteps;

    return states;
}

void
fastpm_tevo_destroy_states(FastPMStates * states)
{
    free(states->table);
    free(states->timesteps);
}

static double
i2t(FastPMStates * states, int i)
{
    int d = i / states->cycle_len;
    double r = (i - states->cycle_len * d) / (1.0 * states->cycle_len);

    if(d >= states->cycles) {
        return states->timesteps[states->cycles];
    } else if(d < 0) {
        return states->timesteps[0];
    }

    if(r != 0.0) {
        return exp((1 - r) * log(states->timesteps[d]) + r * log(states->timesteps[d+1]));
    } else {
        /* no interpolation return an exact value from the table. */
        /* otherwise we mess up comparisons from the intial condition. */
        return states->timesteps[d];
    }
}

void
fastpm_tevo_transition_init(FastPMTransition * transition, FastPMStates * states, int istart, int iend)
{
    FastPMState * start = &states->table[istart];
    FastPMState * end = &states->table[iend];

    transition->states = states;
    transition->istart = istart;
    transition->iend = iend;
    transition->start = start;
    transition->end = end;

    if(start->force != end->force) {
        /* Force */
        transition->action = FASTPM_ACTION_FORCE;
        if(start->x != end->x) {
            fastpm_raise(-1, "A force action must have identical x stamp\n");
        }
        transition->a.i = i2t(states, start->force);
        transition->a.f = i2t(states, end->force);
        transition->a.r = i2t(states, end->x);
        transition->i.i = start->force;
        transition->i.f = end->force;
        transition->i.r = end->x;
    }
    if(start->v != end->v) {
        /* Kick */
        transition->action = FASTPM_ACTION_KICK;
        if(start->force != end->force) {
            fastpm_raise(-1, "A kick action must have identical a stamp\n");
        }
        transition->a.i = i2t(states, start->v);
        transition->a.f = i2t(states, end->v);
        transition->a.r = i2t(states, end->force);
        transition->i.i = start->v;
        transition->i.f = end->v;
        transition->i.r = end->force;
    }
    if(start->x != end->x) {
        /* Drift */
        transition->action = FASTPM_ACTION_DRIFT;
        if(start->v != end->v) {
            fastpm_raise(-1, "A drift action must have identical v stamp\n");
        }
        transition->a.i = i2t(states, start->x);
        transition->a.f = i2t(states, end->x);
        transition->a.r = i2t(states, end->v);
        transition->i.i = start->x;
        transition->i.f = end->x;
        transition->i.r = end->v;
    }
}

int
fastpm_tevo_transition_find_dual(FastPMTransition * transition, FastPMTransition * dual)
{
    /* Find the dual action that updates the dual to the last state, invert it */
    if(transition->end->x != transition->end->v) {
            fastpm_raise(-1, "Only transitions towards a synced x and v has a dual.\n");
    }
    int i;
    FastPMStates * states = transition->states;
    enum FastPMAction dual_action = FASTPM_ACTION_FORCE;
    switch(transition->action) {
        case FASTPM_ACTION_DRIFT:
            dual_action = FASTPM_ACTION_KICK;
        break;
        case FASTPM_ACTION_KICK:
            dual_action = FASTPM_ACTION_DRIFT;
        break;
        default:
            fastpm_raise(-1, "Only Kick and Drift has dual transitions\n");
    }
    for(i = transition->istart; i >= 0; i --) {
        fastpm_tevo_transition_init(dual, states, i - 1, i);
        if(dual->action == dual_action) break;
    }
    if(i == -1) { /* not found */
        return 0;
        fastpm_raise(-1, "Dual transition not found. The state table is likely run. Look at states->table.\n");
    }

    /* the reference is in the future. Need to revert this change */
    fastpm_tevo_transition_init(dual, states, i, i - 1);

    if(dual->a.r != transition->a.i) {
        fastpm_raise(-1, "dual transition reference is not the same as my initial state.\n");
    }
    return 1;
}

int
fastpm_tevo_transition_find_next(FastPMTransition * transition, FastPMTransition * next)
{
    int i;
    FastPMStates * states = transition->states;

    for(i = transition->iend; states->table[i + 1].force != -1; i++) {
        fastpm_tevo_transition_init(next, states, i, i + 1);
        if(next->action == transition->action) return 1;
    }
    return 0;
}

