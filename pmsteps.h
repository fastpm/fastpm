#ifndef STEPPING_H
#define STEPPING_H

void 
stepping_kick(PMStore * pi, PMStore * po,
        double ai, double af, double ac,
        double OmegaM);
void 
stepping_drift(PMStore * pi, PMStore * po,
        double ai, double af, double ac,
        double OmegaM);

void stepping_set_boost(double boost);

int stepping_get_nsteps();

void 
stepping_get_times(int istep,
    double * a_x, double * a_x1,
    double * a_v, double * a_v1);

void 
stepping_set_snapshot(PMStore * pi, PMStore * po,
                double aout, double a_x, double a_v, double OmegaM);

void stepping_init(Parameters * param);
#endif
