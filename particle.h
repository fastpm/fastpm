#ifndef PARTICLE_H
#define PARTICLE_H
#include <stddef.h>
#include <stdint.h>

typedef struct {
    float (* x)[3]; 
    float (* v)[3];  
    int64_t (* id); 
    float (* dx1) [3]; // ZA displacement
    float (* dx2) [3]; // 2LPT displacement

    float (* force)[3]; 
    int np_local;
    int np_allocated;
    size_t np_total;
    float np_average;

/* for a snapshot */
    float a; //, a_x, a_v;
    float boxsize;
    float qfactor;
    int nc;
    float omega_m, h;
    int seed;
    //char filename[64];
    char* filename;
} Particles;

#endif
