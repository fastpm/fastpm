#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
#include <fastpm/timemachine.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"
#include "solver-pm-internal.h"

int fastpm_tevo_block_len(FastPMTEEntry *template) {
    /* Return number of elements in a template */
    int i = 0;
    while(template[i].a != -1) { i++; }
    return i++;
}

FastPMTEStates *fastpm_tevo_generate_states(FastPMTEStates *states, int cycles, FastPMTEEntry *template, double *ts) {
    /* Generate state table */
    int i, j, len = fastpm_tevo_block_len(template), N = len * cycles;
    
    double *timesteps = malloc((cycles + 1) * sizeof(double));
    FastPMTEEntry *table = malloc((N + 3) * sizeof(FastPMTEEntry));
    
    table[0].x = 0; // Initial conditions
    table[0].v = 0;
    table[0].a = -2;
    
    table[1].x = 0; // Force calculation
    table[1].v = 0;
    table[1].a = 0;
    
    for(i = 0; i < cycles; i++) { // Templates
        for(j = 0; j < len; j++) {
            table[j + i * len + 2].a = table[i * len + 1].a + template[j].a;
            table[j + i * len + 2].x = table[i * len + 1].x + template[j].x;
            table[j + i * len + 2].v = table[i * len + 1].v + template[j].v;
        }
    }
    
    table[N+2].a = -1; // End of table
    table[N+2].x = -1;
    table[N+2].v = -1;
    
    states->table = table; // Pack in struct
    states->cycle_len = template[len-1].a;
    states->cycles = cycles;
    
    for(i = 0; i <= cycles; i++) { timesteps[i] = ts[i]; }
    states->timesteps = timesteps;
    
    return states;
}

void fastpm_tevo_destroy_states(FastPMTEStates *states) { free(states->table); free(states->timesteps); }

double fastpm_tevo_i2t(FastPMTEStates *states, int i) {
    int d = i / states->cycle_len;
    double r = (i - states->cycle_len * d) / (1.0 * states->cycle_len);

    if(d >= states->cycles) {
        return states->timesteps[states->cycles];
    } else if(d < 0) {
        return states->timesteps[0];
    }

    return exp((1 - r) * log(states->timesteps[d]) + r * log(states->timesteps[d+1]));
}

void fastpm_tevo_print_states(FastPMTEStates *states) {

    int i = 0;

    printf("# | a | x | v | OP\n");
    printf("%d | \\ | %d | %d | INIT. COND.\n", i, states->table[i].x, states->table[i].v);

    FastPMTEStep lastkick = {-1, -1, -1}, thiskick;
    FastPMTEStep lastdrift = {-1, -1, -1}, thisdrift;

    enum {
        FORCE = 1,
        KICK = 2,
        DRIFT = 3,
    };

    int action = 0;
    i = 1;
    
    while(states->table[i].a != -1) {
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

        printf("%d | %d | %d | %d | ", i, states->table[i].a, states->table[i].x, states->table[i].v);

        if(states->table[i].x == states->table[i].v && states->table[i].x != states->table[i].a) {
            /* Interpolation */
            if(states->table[i].v != states->table[i-1].v) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    thiskick.f, thiskick.i, thiskick.r,
                    lastdrift.i, lastdrift.f, lastdrift.r
                );
            }
            if(states->table[i].x != states->table[i-1].x) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    lastkick.i, lastkick.f, lastkick.r,
                    thisdrift.f, thisdrift.i, thisdrift.r
                );
            }
        }

        switch(action) {
            case FORCE:
                printf("FORCE");
                break;
            case KICK:

                printf("KICK (v=%g, %g, a=%g)", fastpm_tevo_i2t(states, thiskick.f),
                        fastpm_tevo_i2t(states, thiskick.i),
                        fastpm_tevo_i2t(states, thiskick.r));
                lastkick = thiskick;
                break;
            case DRIFT:

                printf("DRIFT(x=%g, %g, v=%g)", fastpm_tevo_i2t(states, thisdrift.f),
                        fastpm_tevo_i2t(states, thisdrift.i),
                        fastpm_tevo_i2t(states, thisdrift.r));
                lastdrift = thisdrift;
                break;
        }

        printf("\n");

        i++;
    }
}

void fastpm_tevo_evolve(FastPMSolver * fastpm, double * time_step, int nstep) 
{
    CLOCK(decompose);
    CLOCK(force);
    CLOCK(kick);
    CLOCK(drift);
    CLOCK(afterforce);
    CLOCK(beforekick);
    CLOCK(beforedrift);
    CLOCK(correction);

    FastPMStore * p = fastpm->p;
    MPI_Comm comm = fastpm->comm;

    pm_2lpt_evolve(time_step[0], p, fastpm->omega_m, fastpm->USE_DX1_ONLY);

    MPI_Barrier(comm);

    CLOCK(warmup);

    if(fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
        /* If doing COLA, v_res = 0 at initial. */
        memset(p->v, 0, sizeof(p->v[0]) * p->np);
    }

    LEAVE(warmup);

    MPI_Barrier(comm);

    /* The last step is the 'terminal' step */
    
    FastPMTEEntry template[] = {
    {0, 0, 1}, /* Kick */
    {0, 1, 1}, /* Drift */
    {0, 2, 1}, /* Drift */
    {2, 2, 1}, /* Force */
    {2, 2, 2}, /* Kick */
    {-1, -1, -1} /* End of table */
    };

    FastPMTEStates *states = malloc(sizeof(FastPMTEStates));
    fastpm_tevo_generate_states(states, nstep-1, template, time_step);
    
    int i = 0;

    FastPMTEStep lastkick = {-1, -1, -1}, thiskick;
    FastPMTEStep lastdrift = {-1, -1, -1}, thisdrift;

    enum {
        FORCE = 1,
        KICK = 2,
        DRIFT = 3,
    };

    int action = 0;
    i = 1;
    
    printf("# | a | x | v | p_x | p_ v | OP\n");
    printf("%d | \\ | %d | %d | %g | %g | INIT. COND.\n", i, states->table[i].x, states->table[i].v, p->a_x, p->a_v);

    while(states->table[i].a != -1) {

        fastpm->info.istep = i;
        fastpm->info.a_x = states->table[i-1].x;
        fastpm->info.a_x1 = states->table[i].x;
        fastpm->info.a_v = states->table[i-1].v;
        fastpm->info.a_v1 = states->table[i].v;
    
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

        printf("%d | %d | %d | %d | %g | %g | ", i, states->table[i].a, states->table[i].x, states->table[i].v, p->a_x, p->a_v);

        if(states->table[i].x == states->table[i].v && states->table[i].x != states->table[i].a) {
            /* Interpolation */
            if(states->table[i].v != states->table[i-1].v) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    thiskick.f, thiskick.i, thiskick.r,
                    lastdrift.i, lastdrift.f, lastdrift.r
                );
            }
            if(states->table[i].x != states->table[i-1].x) {
                printf("I with K(%d %d %d) D(%d %d %d) ",
                    lastkick.i, lastkick.f, lastkick.r,
                    thisdrift.f, thisdrift.i, thisdrift.r
                );
            }
        }

        switch(action) {
            case FORCE:
                printf("FORCE");

                /* Calculate forces */
                VPM * vpm = vpm_find(fastpm->vpm_list, fastpm_tevo_i2t(states, states->table[i].a));
                fastpm->pm = &vpm->pm;
                PM * pm = fastpm->pm;

                fastpm->info.Nmesh = fastpm->pm->init.Nmesh;

                fastpm_painter_init(fastpm->painter, pm, fastpm->PAINTER_TYPE, fastpm->painter_support);

                FastPMFloat * delta_k = pm_alloc(pm);

                ENTER(force);
                fastpm_calculate_forces(fastpm, delta_k);
                LEAVE(force);

                pm_free(pm, delta_k);

                break;
            case KICK:
                printf("KICK (v=%g, %g, a=%g) [v=%d, %d, a=%d]", fastpm_tevo_i2t(states, thiskick.f),
                                                fastpm_tevo_i2t(states, thiskick.i),
                                                fastpm_tevo_i2t(states, thiskick.r),
                                                thiskick.f,
                                                thiskick.i,
                                                thiskick.r);
                lastkick = thiskick;
                /* Do kick */
                ENTER(kick);
                fastpm_kick_store(fastpm, p, p, fastpm_tevo_i2t(states, thiskick.f));
                LEAVE(kick);
                break;
            case DRIFT:
                printf("DRIFT(x=%g, %g, v=%g) [x=%d, %d, v=%d]", fastpm_tevo_i2t(states, thisdrift.f),
                                                fastpm_tevo_i2t(states, thisdrift.i),
                                                fastpm_tevo_i2t(states, thisdrift.r),
                                                thisdrift.f,
                                                thisdrift.i,
                                                thisdrift.r);
                lastdrift = thisdrift;

                /* Do drift */
                ENTER(drift);
                fastpm_drift_store(fastpm, p, p, fastpm_tevo_i2t(states, thisdrift.f));
                LEAVE(drift);

                break;
        }

        printf("\n");

        i++;
    }

    printf("%d | %d | %d | %d | %g | %g\n", i, states->table[i].a, states->table[i].x, states->table[i].v, p->a_x, p->a_v);
}
