FASTPM_BEGIN_DECLS

void
fastpm_apply_smoothing_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double sml);

void
fastpm_apply_lowpass_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double kth);

void
fastpm_apply_decic_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to);

void
fastpm_apply_diff_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) ;

void
fastpm_apply_multiply_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double value);

void
fastpm_apply_laplace_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to);

void
fastpm_apply_any_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, fastpm_fkfunc func, void * data);

void
fastpm_apply_normalize_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to);

void
fastpm_apply_modify_mode_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, ptrdiff_t * mode, double value);

FASTPM_END_DECLS
