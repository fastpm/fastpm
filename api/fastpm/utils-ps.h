FASTPM_BEGIN_DECLS

typedef struct {
    size_t size;
    double *k;
    double *p;

    PM * pm;
    double k0;
    double Volume;
    double *Nmodes;
} FastPMPowerSpectrum;

void
fastpm_powerspectrum_init(FastPMPowerSpectrum * ps, PM * pm, FastPMFloat * delta1_k, FastPMFloat * delta2_k);

void
fastpm_powerspectrum_destroy(FastPMPowerSpectrum * ps);

void
fastpm_powerspectrum_write(FastPMPowerSpectrum * ps, char * filename, double N);

double
fastpm_powerspectrum_large_scale(FastPMPowerSpectrum * ps, int Nmax);

FASTPM_END_DECLS
