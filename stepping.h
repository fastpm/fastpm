#ifndef STEPPING_H
#define STEPPING_H

void stepping_kick(Particles* particles,
        double ai, double af, double ac);
void stepping_drift(Particles* particles, 
        double ai, double af, double ac);

void stepping_set_initial(double aout, Particles * particles);

void stepping_set_snapshot(double aout, double a_x, double a_v, Particles * particles, Particles* snapshot);

void stepping_init(Parameters * param);
#endif
