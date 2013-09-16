#ifndef COLA_H
#define COLA_H 1

//void cola_evolve(Particles* const particles, const float Omega_m, 
//		 const float a0, const float a1, const float a2);

void cola_kick(Particles* const particles, const float Omega_m,
	       const float avel1);
void cola_drift(Particles* const particles, const float Omega_m,
		const float apos1);

void cola_evolve(Particles* const particles, const float Omega_m, 
		 const float apos0, const float apos1,
		 const float avel0, const float avel1);

void set_noncola_initial(Particles const * const particles, Snapshot* const snapshot);

void cola_set_snapshot(const double aout, Particles const * const particles, Snapshot* const snapshot);

#endif
