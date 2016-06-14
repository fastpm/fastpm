FASTPM_BEGIN_DECLS

typedef enum {
    FASTPM_FNL_NONE,
    FASTPM_FNL_LOCAL,
} FastPMPNGaussianType;

typedef struct {
    /* Linear theory power spectrum */
    fastpm_fkfunc pkfunc;
    void * pkdata;

    FastPMPNGaussianType type;
    /* fNL */
    double fNL;
    double kmax_primordial;
    double h;
    double scalar_amp;
    double scalar_spectral_index;
    double scalar_pivot;
    /* FIXME: other parameters that are initialized by the struct. */

    /* private: */
    double Volume;
} FastPMPNGaussian;

void
fastpm_png_induce_correlation(FastPMPNGaussian * png, PM * pm, FastPMFloat * delta_k);

FASTPM_END_DECLS
