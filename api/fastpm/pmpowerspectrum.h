FASTPM_BEGIN_DECLS

typedef struct {
    size_t size;
    double *k;
    double *p;
    double *N;
} PowerSpectrum;

void 
pm_calculate_powerspectrum(PM * pm, FastPMFloat * delta_k, PowerSpectrum * ps);

void 
power_spectrum_init(PowerSpectrum * ps, size_t size);

void 
power_spectrum_destroy(PowerSpectrum * ps);

void
power_spectrum_write(PowerSpectrum * ps, PM * pm, double ntotal, char * basename, int random_seed, double aout);

FASTPM_END_DECLS
