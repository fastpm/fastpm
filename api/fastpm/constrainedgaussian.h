FASTPM_BEGIN_DECLS

typedef struct {
    double x[3];
    double c;
} FastPMConstraint;

typedef struct {
    FastPMConstraint * constraints;
} FastPMConstrainedGaussian;

typedef struct {
    int size;
    double * xi;
    double step_size;
} FastPM2PCF;

double
fastpm_2pcf_eval(FastPM2PCF * self, double r);

void
fastpm_2pcf_from_powerspectrum(FastPM2PCF * self, fastpm_fkfunc pkfunc, void * data, double r_max, int steps);

void
fastpm_cg_apply_constraints(FastPMConstrainedGaussian * cg, PM * pm, FastPM2PCF * xi, FastPMFloat * delta_k);

FASTPM_END_DECLS
