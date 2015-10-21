#ifndef STEPPING_H
#define STEPPING_H
typedef struct {
    double omega_m;
    int mode;
    int stdda;
} PMStepper;

void 
stepping_kick(PMStepper * stepper, 
        PMStore * pi, PMStore * po,
        double ai, double af, double ac);
void 
stepping_drift(PMStepper * stepper, 
        PMStore * pi, PMStore * po,
        double ai, double af, double ac);


void stepping_set_boost(double boost);

void 
stepping_set_snapshot(PMStepper * stepper, 
                PMStore * pi, PMStore * po,
                double aout, double a_x, double a_v);


void stepping_init(PMStepper * stepper, double omega_m, int force_mode, int stdDA);
#endif
