#ifndef STEPPING_H
#define STEPPING_H

void stepping_kick(PMStore* particles,
        double ai, double af, double ac);
void stepping_drift(PMStore* particles, 
        double ai, double af, double ac);
void stepping_set_boost(double boost);

int stepping_get_nsteps();
void stepping_get_times(int istep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1);

void stepping_set_initial(double aout, PMStore * particles, double shift[3]);

void stepping_set_snapshot(double aout, double a_x, double a_v, PMStore * particles, PMStore* snapshot);

void stepping_init(Parameters * param);
#endif
