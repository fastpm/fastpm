#ifndef COLA_H
#define COLA_H 1

void cola_kick(Particles* const particles, const float Omega_m,
        const float ai, const float af, const float ac);
void cola_drift(Particles* const particles, const float Omega_m,
        const float ai, const float af, const float ac);

void set_noncola_initial(Particles const * const particles, Snapshot* const snapshot);

void cola_set_snapshot(const double aout, Particles const * const particles, Snapshot* const snapshot);

#endif
