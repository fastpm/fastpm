#ifndef PM_H
#define PM_H 1

#include "particle.h"

void pm_set_diff_order(int order);
void pm_init(double boxsize, int nc);
void pm_set_size(int nc_pm_factor);

void pm_free(void);

void pm_calculate_forces(Particles*);
double * pm_compute_power_spectrum(size_t * nk);

#endif
