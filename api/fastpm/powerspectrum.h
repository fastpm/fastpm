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

int
fastpm_powerspectrum_init_from_string(FastPMPowerSpectrum * ps, const char * string);

void
fastpm_powerspectrum_destroy(FastPMPowerSpectrum * ps);

void
fastpm_powerspectrum_write(FastPMPowerSpectrum * ps, char * filename, double N);

double
fastpm_powerspectrum_large_scale(FastPMPowerSpectrum * ps, int Nmax);

double
fastpm_powerspectrum_eval(FastPMPowerSpectrum * ps, double k);

double
fastpm_powerspectrum_sigma(FastPMPowerSpectrum * ps, double R);

void
fastpm_powerspectrum_scale(FastPMPowerSpectrum * ps, double factor);

FASTPM_END_DECLS
