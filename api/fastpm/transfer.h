void 
fastpm_apply_smoothing_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double sml);

void 
fastpm_apply_diff_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) ;

void fastpm_apply_za_hmc_force_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir);
