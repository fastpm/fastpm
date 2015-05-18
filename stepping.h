#ifndef STEPPING_H
#define STEPPING_H

void stepping_kick(Particles* particles, double Omega_m,
        double ai, double af, double ac);
void stepping_drift(Particles* particles, double Omega_m,
        double ai, double af, double ac);

void stepping_set_subtract_lpt(int flag);
void stepping_set_std_da(int flag);
void stepping_set_no_pm(int flag);
void stepping_set_initial(double aout, double omega_m, Particles * particles);

void stepping_set_snapshot(double aout, double a_x, double a_v, Particles * particles, Particles* snapshot);

#endif
