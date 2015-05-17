#ifndef STEPPING_H
#define STEPPING_H

void stepping_kick(Particles* const particles, const double Omega_m,
        const double ai, const double af, const double ac);
void stepping_drift(Particles* const particles, const double Omega_m,
        const double ai, const double af, const double ac);

void stepping_set_subtract_lpt(int flag);
void stepping_set_std_da(int flag);
void stepping_set_no_pm(int flag);
void stepping_set_initial(const double aout, const double omega_m, Particles const * const particles);

void stepping_set_snapshot(const double aout, double a_x, double a_v, Particles const * const particles, Snapshot* const snapshot);

#endif
