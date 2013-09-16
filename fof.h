#ifndef FOF_H
#define FOF_H 1

#include "particle.h"

size_t fof_calc_memory(const int np_alloc, const int nc);
void fof_init(const int np_alloc, const int nc, void* mem, size_t mem_size);
void fof_find_halos(Snapshot* snapshot, const float ll);
//void fof_write_halos(char filename[], void* buf, size_t size);
void fof_write_halos(char filename[]);

#endif
