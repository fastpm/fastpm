#ifndef PM_H
#define PM_H 1

#include "particle.h"
typedef double (*pm_mond_func)(double k2);
double PM_NO_MOND(double k2);

void pm_set_diff_order(int order);
void pm_set_mond(pm_mond_func mond_mu_k);

void pm_init(const int nc_pm, const int nc_pm_factor, const float boxsize,
         int many);
//void move_particles2(Particles*);
void pm_finalize(void);
void pm_calculate_forces(Particles*);
double * pm_compute_power_spectrum(size_t * nk);

//void PtoMesh(const Particle Pz[], const int NumPart);
//void forces();
//void pm_exchange_particles(Particles* const particles);

//int pm_set_cic_density(Particle* const p, int np_local, const int np_alloc);

#endif
