#ifndef PARTICLE_H
#define PARTICLE_H

typedef float float3[3];

typedef struct {
  long long id;
  float x[3];
  float v[3];   // velocity
  float dx1[3]; // ZA displacement
  float dx2[3]; // 2LPT displacement
} Particle;

typedef struct {
  Particle* p;
  float3* force;

  int np_local, np_allocated;
  long long np_total;
  float np_average;
} Particles;

typedef struct {
  long long id;
  float x[3];
  float v[3];
  int iGroup;
  int OriginalTask;
  long long OriginalIndex;
} ParticleMinimum;

typedef struct {
  ParticleMinimum* p;
  int np_local;
  int np_allocated;
  long long np_total;
  float np_average;
  float a; //, a_x, a_v;
  float boxsize;
  float qfactor;
  int nc;
  float omega_m, h;
  int seed;
  //char filename[64];
  char* filename;
} Snapshot;

#endif
