FASTPM_BEGIN_DECLS

typedef struct {
    /* in units of real numbers, not bytes. */
    ptrdiff_t start[3];
    ptrdiff_t size[3];
    ptrdiff_t strides[3];
    ptrdiff_t total;
} PMRegion;

/* 
 * Allocate memory for FFT/painting in PM. 
 * */

FastPMFloat * pm_alloc(PM * pm);
void pm_free(PM * pm, FastPMFloat * buf);
void pm_assign(PM * pm, FastPMFloat * from, FastPMFloat * to);

/* property accessors of PM objects */
size_t pm_allocsize(PM * pm);
MPI_Comm pm_comm(PM * pm);
double pm_norm(PM * pm);
ptrdiff_t * pm_nmesh(PM * pm);
int * pm_nproc(PM * pm);
double * pm_boxsize(PM * pm);
double pm_volume(PM * pm);
PMRegion * pm_i_region(PM * pm);
PMRegion * pm_o_region(PM * pm);

void pm_unravel_o_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]);
void pm_unravel_i_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]);

ptrdiff_t pm_ravel_o_index(PM * pm, ptrdiff_t i[3]);
ptrdiff_t pm_ravel_i_index(PM * pm, ptrdiff_t i[3]);

int pm_inc_o_index(PM * pm, ptrdiff_t i[3]);
int pm_inc_i_index(PM * pm, ptrdiff_t i[3]);


typedef struct {
    ptrdiff_t start;
    ptrdiff_t end;
    ptrdiff_t ind;
    ptrdiff_t i[3];
    ptrdiff_t iabs[3];

    float *k_finite[3]; /* k, 4 point central */
    float *k[3]; /* k */
    float *kk_finite[3]; /* k ** 2, 3 point central */
    float *kk_finite2[3]; /* k ** 2, 5 point central */
    float *kk[3];  /* k ** 2 */

    PM * pm;
} PMKIter;

void 
pm_kiter_init(PM * pm, PMKIter * iter);

int pm_kiter_stop(PMKIter * iter);

void pm_kiter_next(PMKIter * iter);

double pm_kiter_get_kmag(PMKIter * iter);

typedef struct {
    ptrdiff_t start;
    ptrdiff_t end;
    ptrdiff_t ind;
    ptrdiff_t i[3];
    ptrdiff_t iabs[3];

    PM * pm;
} PMXIter;

void 
pm_xiter_init(PM * pm, PMXIter * iter);

int pm_xiter_stop(PMXIter * iter);

void pm_xiter_next(PMXIter * iter);

/* 
 * r2c is out-of-place and c2r is in-place.
 * */
void 
pm_r2c(PM * pm, FastPMFloat * from, FastPMFloat * to);
void 
pm_c2r(PM * pm, FastPMFloat * inplace);


FASTPM_END_DECLS
