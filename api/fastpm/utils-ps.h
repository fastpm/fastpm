FASTPM_BEGIN_DECLS

typedef struct {
    PM * pm;
    size_t size;
    double Volume;
    double *k;
    double *p;
    double *Nmodes;
} FastPMPowerSpectrum;

void
fastpm_powerspectrum_init(FastPMPowerSpectrum * ps, PM * pm);

void
fastpm_powerspectrum_measure(FastPMPowerSpectrum * ps, FastPMFloat * delta1_k, FastPMFloat * delta2_k);

void
fastpm_powerspectrum_destroy(FastPMPowerSpectrum * ps);

void
fastpm_powerspectrum_write(FastPMPowerSpectrum * ps, char * filename, double N);

double
fastpm_calculate_large_scale_power(PM * pm, FastPMFloat * delta_k, int Nmax);

FASTPM_END_DECLS
