FASTPM_BEGIN_DECLS

/* 
 * Allocate memory for FFT/painting in PM. 
 * */

FastPMFloat * pm_alloc(PM * pm);
void pm_free(PM * pm, FastPMFloat * buf);
void pm_assign(PM * pm, FastPMFloat * from, FastPMFloat * to);

/* property accessors of PM objects */
size_t pm_size(PM * pm);
ptrdiff_t * pm_nmesh(PM * pm);
double * pm_boxsize(PM * pm);
PMRegion * pm_i_region(PM * pm);
PMRegion * pm_o_region(PM * pm);

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */

void
fastpm_utils_paint(PM * pm, PMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k);

void 
fastpm_utils_dump(PM * pm , char * filename, FastPMFloat *data);

FASTPM_END_DECLS

