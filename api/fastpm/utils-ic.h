typedef double (*fastpm_pkfunc)(double k, void * data);

enum FastPMFillDeltaKScheme {
    FASTPM_DELTAK_GADGET,
    FASTPM_DELTAK_FAST,
    FASTPM_DELTAK_SLOW,
};
void 
fastpm_utils_fill_deltak(PM * pm, FastPMFloat * delta_k, 
        int seed, fastpm_pkfunc pk, void * pkdata, enum FastPMFillDeltaKScheme scheme);

void 
fastpm_utils_induce_correlation(PM * pm, FastPMFloat * g_x, FastPMFloat * delta_k, fastpm_pkfunc pk, void * pkdata);

