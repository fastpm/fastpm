#include <string.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>
#include "pmpfft.h"

FastPMFloat * pm_alloc_details(PM * pm, const char * file, const int line)
{
    void * p = fastpm_memory_alloc_details(pm->mem, "PMAlloc", sizeof(FastPMFloat) * pm->allocsize, FASTPM_MEMORY_HEAP, file, line);
    memset(p, 0, pm->allocsize * sizeof(FastPMFloat));
    return p;
}

void 
pm_free(PM * pm, FastPMFloat * data)
{
    fastpm_memory_free(pm->mem, data);
}

void 
pm_assign(PM * pm, FastPMFloat * from, FastPMFloat * to) 
{
    memcpy(to, from, sizeof(from[0]) * pm->allocsize);
}

void 
pm_clear(PM * pm, FastPMFloat * buf)
{
    memset(buf, 0, sizeof(buf[0]) * pm->allocsize);
}

size_t 
pm_allocsize(PM * pm)
{
    return pm->allocsize;
}

double 
pm_norm(PM * pm)
{
    return pm->Norm;
}

MPI_Comm pm_comm(PM * pm)
{
    return pm->Comm2D;
}

ptrdiff_t * pm_nmesh(PM * pm) {
    return pm->Nmesh;
}

int * pm_nproc(PM * pm) {
    return pm->Nproc;
}

double * pm_boxsize(PM * pm) {
    return pm->BoxSize;
}

double pm_volume(PM * pm) {
    return pm->Volume;
}

int
pm_unbalanced(PM * pm)
{
    if(pm->Nmesh[0] % pm->Nproc[0] != 0 ||
       pm->Nmesh[1] % pm->Nproc[1] != 0) {
        return 1;
    }
    return 0;
}


PMRegion * pm_i_region(PM * pm) {
    return &pm->IRegion;
}

PMRegion * pm_o_region(PM * pm) {
    return &pm->ORegion;
}

static void 
pm_create_k_factors(PM * pm, PMKIter * iter);

static void 
pm_destroy_k_factors(PMKIter * iter);

void 
pm_kiter_init(PM * pm, PMKIter * iter) 
{
    /* static schedule the openmp loops. start, end is in units of 'real' numbers.
     *
     * i is in units of complex numbers.
     *
     * We call pm_unravel_o_index to set the initial i[] for each threads,
     * then rely on pm_inc_o_index to increment i, because the former is 
     * much slower than pm_inc_o_index and would eliminate threading advantage.
     *
     * */
#ifdef _OPENMP
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();
#else
    int nth = 1;
    int ith = 0;
#endif
    iter->start = ith * pm->ORegion.total / nth * 2;
    iter->end = (ith + 1) * pm->ORegion.total / nth * 2;

    /* do not unravel if we are not looping at all. 
     * This fixes a FPE when
     * the rank has ORegion.total == 0 
     * -- with PFFT the last transposed dimension
     * on some ranks will be 0 */
    if(iter->end > iter->start) 
        pm_unravel_o_index(pm, iter->start / 2, iter->i);

    iter->ind = iter->start;
    iter->pm = pm;

    int d;
    for(d = 0; d < 3; d ++) {
        iter->iabs[d] = iter->i[d] + pm->ORegion.start[d];
    }

    pm_create_k_factors(pm, iter);

}

double pm_kiter_get_kmag(PMKIter * iter)
{
    int d = 0;
    double kk = 0;
    for(d = 0; d < 3; d ++) {
        kk += iter->kk[d][iter->iabs[d]];
    }
    return sqrt(kk);
}

int pm_kiter_stop(PMKIter * iter) 
{
    int stop = !(iter->ind < iter->end);
    if(stop) {
        pm_destroy_k_factors(iter);
    }
    return stop;
}

void pm_kiter_next(PMKIter * iter) 
{
    iter->ind += 2;
    pm_inc_o_index(iter->pm, iter->i);
    int d;
    for(d = 0; d < 3; d ++) {
        iter->iabs[d] = iter->i[d] + iter->pm->ORegion.start[d];
    }
}

void 
pm_xiter_init(PM * pm, PMXIter * iter) 
{
#ifdef _OPENMP
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();
#else
    int nth = 1;
    int ith = 0;
#endif
    iter->start = ith * pm->IRegion.total / nth;
    iter->end = (ith + 1) * pm->IRegion.total / nth;

    /* do not unravel if we are not looping at all. 
     * This fixes a FPE when
     * the rank has ORegion.total == 0 
     * -- with PFFT the last transposed dimension
     * on some ranks will be 0 */
    if(iter->end > iter->start) 
        pm_unravel_i_index(pm, iter->start, iter->i);

    iter->ind = iter->start;
    iter->pm = pm;
    int d;
    for(d = 0; d < 3; d ++) {
        iter->iabs[d] = iter->i[d] + pm->IRegion.start[d];
    }
}

int pm_xiter_stop(PMXIter * iter) 
{
    int stop = !(iter->ind < iter->end);

    return stop;
}

void pm_xiter_next(PMXIter * iter) 
{
    iter->ind += pm_inc_i_index(iter->pm, iter->i);
    int d;
    for(d = 0; d < 3; d ++) {
        iter->iabs[d] = iter->i[d] + iter->pm->IRegion.start[d];
    }
}

static double 
sinc_unnormed(double x);

static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}

static double 
diff_kernel(double w) 
{
    /* order N = 1 super lanzcos kernel */
    /* 
     * This is the same as GADGET-2 but in fourier space: 
     * see gadget-2 paper and Hamming's book.
     * c1 = 2 / 3, c2 = 1 / 12
     * */
    return 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}

static void 
pm_create_k_factors(PM * pm, PMKIter * iter) 
{ 
    /* This function populates fac with precalculated values that
     * are useful for force calculation. 
     * e.g. k**2 and the finite differentiation kernels. 
     * precalculating them means in the true kernel we only need a 
     * table look up. watch out for the offset ORegion.start
     * */
    int d;
    ptrdiff_t ind;
    for(d = 0; d < 3; d++) {
        double CellSize = pm->BoxSize[d] / pm->Nmesh[d];

        iter->k[d] = malloc(sizeof(float) * pm->Nmesh[d]);
        iter->k_finite[d] = malloc(sizeof(float) * pm->Nmesh[d]);
        iter->kk[d] = malloc(sizeof(float) * pm->Nmesh[d]);
        iter->kk_finite[d] = malloc(sizeof(float) * pm->Nmesh[d]);
        iter->kk_finite2[d] = malloc(sizeof(float) * pm->Nmesh[d]);

        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            float k = pm->MeshtoK[d][ind];
            float w = k * CellSize;
            float ff1 = sinc_unnormed(0.5 * w);
            float ff2 = sinc_unnormed(w);

            iter->k[d][ind] = k;
            iter->kk[d][ind] = k * k;
/* 4 point central diff */
            iter->k_finite[d][ind] = 1 / CellSize * diff_kernel(w);
/* naive */
//            iter->k_finite[d][ind] = k;

/* 5 point central diff */
            iter->kk_finite2[d][ind] = k * k * ( 4 / 3.0 * ff1 * ff1 - 1 / 3.0 * ff2 * ff2);
/* 3 point central diff */
            iter->kk_finite[d][ind] = k * k * (ff1 * ff1);
/* naive */
//            iter->kk_finite[d][ind] = k * k;
        }
    } 
}

void 
pm_destroy_k_factors(PMKIter * iter) 
{
    int d;
    for(d = 0; d < 3; d ++) {
        free(iter->k_finite[d]);
        free(iter->k[d]);
        free(iter->kk_finite2[d]);
        free(iter->kk_finite[d]);
        free(iter->kk[d]);
    }
}

double
pm_compute_variance(PM * pm, FastPMFloat * complx)
{
    PMKIter kiter;
    double v = 0;
    for(pm_kiter_init(pm, &kiter);
        !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        int w = 2;
        if (kiter.iabs[2] == 0 || kiter.iabs[2] == pm->Nmesh[2] / 2) {
            w = 1;
        }
        v += w * complx[kiter.ind + 0] * complx[kiter.ind + 0];
        v += w * complx[kiter.ind + 1] * complx[kiter.ind + 1];
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    v = v / pm_norm(pm);
    return v;
}

PM *
fastpm_create_pm(int Ngrid, int NprocY, int transposed, double BoxSize, MPI_Comm comm)
{
    PM * pm = malloc(sizeof(PM));

    PMInit pminit = {
        .Nmesh = Ngrid,
        .BoxSize = BoxSize,
        .NprocY = NprocY,
        .transposed = transposed,
        .use_fftw = 0,
    };

    pm_init(pm, &pminit, comm);
    return pm;
}


void
fastpm_free_pm(PM * pm)
{
    pm_destroy(pm);
    free(pm);
}

void
pm_check_values(PM * pm, FastPMFloat * field, const char * fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);
    char * buffer = fastpm_strdup_vprintf(fmt, argp);
    va_end(argp);

    ptrdiff_t i;
    size_t oo = 0;
    for(i = 0; i < pm->allocsize; i ++) {
        FastPMFloat value = field[i];
        if (value > 1e15 || value < -1e15 || value != value) {
            oo ++;
        }
    }
    if(oo != 0) {
        fastpm_ilog(INFO, "%s: Task %d has %td field values that are out of bounds\n",
                buffer, pm->ThisTask, oo);
    }
    free(buffer);
}
