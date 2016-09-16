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

int fastpm_tevo_block_len(FastPMTEEntry *template) {
    /* Return number of elements in a template */
    int i = 0;
    while(template[i].a != -1) { i++; }
    return i++;
}

FastPMTEStates *fastpm_tevo_generate_states(FastPMTEStates *states, int cycles, FastPMTEEntry *template, double *ts) {
    /* Generate state table */
    int i, j, len = fastpm_tevo_block_len(template), N = len * cycles;
    
    double *timesteps = malloc(cycles * sizeof(double));
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
    
    for(i = 0; i < cycles; i++) { timesteps[i] = ts[i]; }
    states->timesteps = timesteps;
    
    return states;
}

void fastpm_tevo_destroy_states(FastPMTEStates *states) { free(states->table); free(states->timesteps); }

double fastpm_tevo_i2t(FastPMTEStates *states, int i) {
    int d = i / states->cycle_len;
    double r = (i - states->cycle_len * d) / (1.0 * states->cycle_len);
    
    if(d+1 >= states->cycles) {
		return states->timesteps[states->cycles - 1];
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
