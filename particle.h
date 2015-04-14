#ifndef PARTICLE_H
#define PARTICLE_H

typedef float float3[3];

typedef struct {
  float x[3];
  union {
      float dx1[3]; // ZA displacement
      int OriginalTask;
  };
  union {
      float dx2[3]; // 2LPT displacement
      int OriginalIndex;
  };
  float v[3];   // velocity
  long long id;
} Particle;

typedef struct {
  Particle* p;
  float3* force;

  int np_local, np_allocated;
  long long np_total;
  float np_average;
} Particles;

typedef struct {
  float x[3];
  float v[3];
  long long id;
} ParticleMinimum;

typedef struct {
  ParticleMinimum* p;
  int np_local;
  int np_allocated;
  long long np_total;
  float np_average;
  float a; //, a_x, a_v;
  float boxsize;
  int nc;
  float omega_m, h;
  int seed;
  //char filename[64];
  char* filename;
} Snapshot;

#endif
