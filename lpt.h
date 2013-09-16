#ifndef LPT_H
#define LPT_H

#include "particle.h"

void lpt_init(const int nc, const void* mem, const size_t size);
void lpt_set_displacement(const double InitTime, const double omega_m, const int Seed, const double Box, Particles* const particles);
int lpt_get_local_nx(void);

#endif
