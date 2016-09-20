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

static int
fastpm_tevo_block_len(FastPMTEEntry *template)
{
    /* Return number of elements in a template */
    int i = 0;
    while(template[i].a != -1) { i++; }
    return i++;
}

FastPMTEStates *
fastpm_tevo_generate_states(FastPMTEStates *states, int cycles, FastPMTEEntry *template, double *ts)
{
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

void
fastpm_tevo_destroy_states(FastPMTEStates * states)
{
    free(states->table);
    free(states->timesteps);
}

static double
fastpm_tevo_i2t(FastPMTEStates * states, int i)
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
fastpm_tevo_transition_init(FastPMTETransition * transition, FastPMTEStates * states, int i, int r, int f)
{
    transition->i = i;
    transition->r = r;
    transition->f = f;
    transition->a_i = fastpm_tevo_i2t(states, i);
    transition->a_r = fastpm_tevo_i2t(states, r);
    transition->a_f = fastpm_tevo_i2t(states, f);
}
