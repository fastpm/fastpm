FASTPM_BEGIN_DECLS

typedef struct {
    /* Linear theory power spectrum */
    fastpm_fkfunc pkfunc;
    void * pkdata;

    /* fNL */
    double fNL;
    /* FIXME: other parameters that are initialized by the struct. */

    /* private: */
    double Volume;
} FastPMPNGaussian;

void
fastpm_png_induce_correlation(FastPMPNGaussian * png, PM * pm, FastPMFloat * delta_k);

FASTPM_END_DECLS
