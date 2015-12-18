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

void
fastpm_paint(PM * pm, PMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k);

void 
fastpm_dump(PM * pm , char * filename, FastPMFloat *data);
