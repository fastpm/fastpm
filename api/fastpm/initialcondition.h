FASTPM_BEGIN_DECLS

enum FastPMFillDeltaKScheme {
    FASTPM_DELTAK_GADGET,
    FASTPM_DELTAK_FAST,
    FASTPM_DELTAK_SLOW,
};
void
fastpm_utils_fill_deltak(PM * pm, FastPMFloat * delta_k,
        int seed, enum FastPMFillDeltaKScheme scheme);

void
fastpm_utils_induce_correlation(PM * pm, FastPMFloat * delta_k, fastpm_pkfunc pk, void * pkdata);

void
fastpm_utils_remove_cosmic_variance(PM * pm, FastPMFloat * delta_k);

FASTPM_END_DECLS
