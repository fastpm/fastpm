#ifndef MEM_H
#define MEM_H 1

#include "particle.h"

typedef struct {
  void *mem1, *mem2;
  size_t size1, size2;
} Memory;


Particles* allocate_particles(const int nc, const int nx, double np_alloc_factor);
Snapshot* allocate_snapshot(const int nc, const int nx, const int np_alloc, void* const mem, const size_t mem_size);

void allocate_shared_memory(const int nc, const int nc_factor, const double np_alloc_factor, Memory* const mem);

#endif
