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

typedef struct {
    float k_finite; /* k, finite */
    float k; /* k */
    float kk_finite; /* k ** 2, on a mesh */
    float kk;  /* k ** 2 */
    float cic;  /* 1 - 2 / 3 sin^2 ( 0.5 k L / N)*/
    float extra;  /* any temporary variable that can be useful. */
} PMKFactors;

typedef struct {
    ptrdiff_t start;
    ptrdiff_t end;
    ptrdiff_t ind;
    ptrdiff_t i[3];
    ptrdiff_t iabs[3];
    PMKFactors * fac[3];
    PM * pm;
} PMKIter;

void 
pm_kiter_init(PM * pm, PMKIter * iter);

int pm_kiter_stop(PMKIter * iter);

void pm_kiter_next(PMKIter * iter);

FASTPM_END_DECLS
