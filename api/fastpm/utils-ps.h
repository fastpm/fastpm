FASTPM_BEGIN_DECLS

typedef struct {
    size_t size;
    double ntotal;
    double *k;
    double *p;
    double *N;
} FastPMPowerSpectrum;

void 
fastpm_utils_calculate_powerspectrum(PM * pm, FastPMFloat * delta_k, FastPMPowerSpectrum * ps, double Ntot);

void 
fastpm_power_spectrum_init(FastPMPowerSpectrum * ps, size_t size);

void 
fastpm_power_spectrum_destroy(FastPMPowerSpectrum * ps);

void
fastpm_power_spectrum_write(FastPMPowerSpectrum * ps, PM * pm, char * filename);

double
fastpm_calculate_large_scale_power(PM * pm, FastPMFloat * delta_k, int Nmax);

FASTPM_END_DECLS
