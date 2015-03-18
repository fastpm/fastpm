#ifndef PM_H
#define PM_H 1

#include "particle.h"

void pm_init(const int nc_pm, const int nc_pm_factor, const float boxsize,
	     void* const mem1, const size_t size1,
	     void* const mem2, const size_t size2,
         int many);
//void move_particles2(Particles*);
void pm_calculate_forces(Particles*);
double * pm_compute_power_spectrum(size_t * nk);

//void PtoMesh(const Particle Pz[], const int NumPart);
//void forces();
//void pm_exchange_particles(Particles* const particles);

//int pm_set_cic_density(Particle* const p, int np_local, const int np_alloc);

#endif
