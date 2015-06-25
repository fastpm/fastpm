#ifndef PM_H
#define PM_H 1

#include "particle.h"

typedef double (* pm_modulator)(double k2, void * data);
void pm_set_mond(pm_modulator mond, void * data);

void pm_set_diff_order(int order);
void pm_init(double boxsize, int nc);
void pm_set_size(int nc_pm_factor);

void pm_free(void);

void pm_calculate_forces(Particles*);
double * pm_compute_power_spectrum(size_t * nk);

double pm_get_broadband();
void pm_enforce_broadband(double power);

#endif
