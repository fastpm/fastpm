FASTPM_BEGIN_DECLS

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */

double 
fastpm_utils_powerspec_white(double k, double * amplitude); /* white noise. */

void
fastpm_utils_paint(PM * pm, FastPMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k,
        fastpm_posfunc get_position, int attribute);

void
fastpm_utils_readout(PM * pm, FastPMStore * p, FastPMFloat * delta_x,
        fastpm_posfunc get_position, int attribute);

void
fastpm_utils_smooth(PM * pm, FastPMFloat * delta_x, FastPMFloat * delta_smooth, double sml);

void 
fastpm_utils_dump(PM * pm , const char * filename, FastPMFloat *data);

void 
fastpm_utils_load(PM * pm , const char * filename, FastPMFloat *data);

double 
fastpm_utils_get_random(uint64_t id);

FASTPM_END_DECLS

