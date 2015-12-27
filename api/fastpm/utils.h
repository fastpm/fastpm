FASTPM_BEGIN_DECLS

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */

void
fastpm_utils_paint(PM * pm, PMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k);

void 
fastpm_utils_dump(PM * pm , char * filename, FastPMFloat *data);

FASTPM_END_DECLS

