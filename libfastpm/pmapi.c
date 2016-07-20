#include <string.h>
#include <mpi.h>
#include <omp.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"

FastPMFloat * pm_alloc(PM * pm)
{
    void * p = fastpm_memory_alloc(pm->mem, sizeof(FastPMFloat) * pm->allocsize, FASTPM_MEMORY_HEAP);
    fastpm_memory_tag(pm->mem, p, "FastPMFloat");
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
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();

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
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();

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

