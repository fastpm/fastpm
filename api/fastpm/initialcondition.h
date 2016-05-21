FASTPM_BEGIN_DECLS

enum FastPMFillDeltaKScheme {
    FASTPM_DELTAK_GADGET,
    FASTPM_DELTAK_FAST,
    FASTPM_DELTAK_SLOW,
};
void
fastpm_ic_fill_gaussiank(PM * pm, FastPMFloat * delta_k,
        int seed, enum FastPMFillDeltaKScheme scheme);

void
fastpm_ic_remove_variance(PM * pm, FastPMFloat * delta_k);

void
fastpm_ic_induce_correlation(PM * pm, FastPMFloat * delta_k, fastpm_fkfunc pk, void * pkdata);

FASTPM_END_DECLS
