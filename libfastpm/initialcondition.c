#include <string.h>
#include <math.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"

/* The following functions fill the gaussian field*/
static void
pmic_fill_gaussian_gadget(PM * pm, FastPMFloat * delta_k, int seed);
static void
pmic_fill_gaussian_fast(PM * pm, FastPMFloat * delta_k, int seed);
static void
pmic_fill_gaussian_slow(PM * pm, FastPMFloat * delta_k, int seed);

void
fastpm_ic_fill_gaussiank(PM * pm, FastPMFloat * delta_k, int seed, enum FastPMFillDeltaKScheme scheme)
{

    /* clear the memory to avoid any modes that we forget to set. */
    memset(delta_k, 0, pm_allocsize(pm) * sizeof(delta_k[0]));

    switch(scheme) {
        case FASTPM_DELTAK_GADGET:
            pmic_fill_gaussian_gadget(pm, delta_k, seed);
            break;
        case FASTPM_DELTAK_FAST:
            pmic_fill_gaussian_fast(pm, delta_k, seed);
            break;
        case FASTPM_DELTAK_SLOW:
            pmic_fill_gaussian_slow(pm, delta_k, seed);
            break;
        default:
            pmic_fill_gaussian_gadget(pm, delta_k, seed);
            break;
    }
}

struct PofK {
    fastpm_fkfunc func;
    void * data;
    double Volume;
} ;

static double _powerspec_to_transfer(double k, struct PofK * pk)
{
    double f = sqrt(pk->func(k, pk->data));
    f *= sqrt(1.0 / pk->Volume);
    return f;
}

void
fastpm_ic_induce_correlation(PM * pm, FastPMFloat * delta_k, fastpm_fkfunc pkfunc, void * data)
{
    struct PofK pk;
    pk.func = pkfunc;
    pk.data = data;
    pk.Volume = pm->Volume;

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) _powerspec_to_transfer, &pk);
}

void
fastpm_ic_remove_variance(PM * pm, FastPMFloat * delta_k)
{
#pragma omp parallel
    {
        PMKIter kiter;

        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k2 = 0;
            int d;
            for(d = 0; d < 3; d++) {
                k2 += kiter.kk[d][kiter.iabs[d]];
            }

            /* https://en.wikipedia.org/wiki/Atan2 */
            double a = delta_k[kiter.ind];
            double b = delta_k[kiter.ind + 1];

            if(a == 0 && b == 0)   {
                delta_k[kiter.ind + 0] = 0;
                delta_k[kiter.ind + 1] = 0;
            } else {
                double phase = atan2(b, a);
                delta_k[kiter.ind + 0] = cos(phase);
                delta_k[kiter.ind + 1] = sin(phase);
            }

        }

    }
}


static inline void 
SETSEED(PM * pm, unsigned int * table[2][2], int i, int j, gsl_rng * rng) 
{ 
    unsigned int seed = 0x7fffffff * gsl_rng_uniform(rng); 

    int ii[2] = {i, (pm->Nmesh[0] - i) % pm->Nmesh[0]};
    int jj[2] = {j, (pm->Nmesh[1] - j) % pm->Nmesh[1]};
    int d1, d2;
    for(d1 = 0; d1 < 2; d1++) {
        ii[d1] -= pm->ORegion.start[0];
        jj[d1] -= pm->ORegion.start[1];
    }
    for(d1 = 0; d1 < 2; d1++)
    for(d2 = 0; d2 < 2; d2++) {
        if( ii[d1] >= 0 && 
            ii[d1] < pm->ORegion.size[0] &&
            jj[d2] >= 0 &&
            jj[d2] < pm->ORegion.size[1]
        ) {
            table[d1][d2][ii[d1] * pm->ORegion.size[1] + jj[d2]] = seed;
        }
    }
}
static inline unsigned int 
GETSEED(PM * pm, unsigned int * table[2][2], int i, int j, int d1, int d2) 
{
    i -= pm->ORegion.start[0];
    j -= pm->ORegion.start[1];
    if(i < 0) abort();
    if(j < 0) abort();
    if(i >= pm->ORegion.size[0]) abort();
    if(j >= pm->ORegion.size[1]) abort();
    return table[d1][d2][i * pm->ORegion.size[1] + j];
}

static void 
SAMPLE(gsl_rng * rng, double * ampl, double * phase)
{
    *phase = gsl_rng_uniform(rng) * 2 * M_PI;
    *ampl = 0;
    do *ampl = gsl_rng_uniform(rng); while(*ampl == 0);
}

static void
pmic_fill_gaussian_gadget(PM * pm, FastPMFloat * delta_k, int seed)
{
    /* Fill delta_k with gadget scheme */
    int d;
    int i, j, k;

    memset(delta_k, 0, sizeof(delta_k[0]) * pm->allocsize);

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(rng, seed);

    unsigned int * seedtable[2][2];
    for(i = 0; i < 2; i ++)
    for(j = 0; j < 2; j ++) {
            seedtable[i][j] = calloc(pm->ORegion.size[0] * pm->ORegion.size[1], sizeof(int));
    }

    for(i = 0; i < pm->Nmesh[0] / 2; i++) {
        for(j = 0; j < i; j++) SETSEED(pm, seedtable, i, j, rng);
        for(j = 0; j < i + 1; j++) SETSEED(pm, seedtable, j, i, rng);
        for(j = 0; j < i; j++) SETSEED(pm, seedtable, pm->Nmesh[0] - 1 - i, j, rng);
        for(j = 0; j < i + 1; j++) SETSEED(pm, seedtable, pm->Nmesh[1] - 1 - j, i, rng);
        for(j = 0; j < i; j++) SETSEED(pm, seedtable, i, pm->Nmesh[1] - 1 - j, rng);
        for(j = 0; j < i + 1; j++) SETSEED(pm, seedtable, j, pm->Nmesh[0] - 1 - i, rng);
        for(j = 0; j < i; j++) SETSEED(pm, seedtable, pm->Nmesh[0] - 1 - i, pm->Nmesh[1] - 1 - j, rng);
        for(j = 0; j < i + 1; j++) SETSEED(pm, seedtable, pm->Nmesh[1] - 1 - j, pm->Nmesh[0] - 1 - i, rng);
    }
    gsl_rng_free(rng);

    ptrdiff_t irel[3];
    for(i = pm->ORegion.start[0]; 
        i < pm->ORegion.start[0] + pm->ORegion.size[0]; 
        i ++) {

        gsl_rng * lower_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
        gsl_rng * this_rng = gsl_rng_alloc(gsl_rng_ranlxd1);

        int ci = pm->Nmesh[0] - i;
        if(ci >= pm->Nmesh[0]) ci -= pm->Nmesh[0];

        for(j = pm->ORegion.start[1]; 
            j < pm->ORegion.start[1] + pm->ORegion.size[1]; 
            j ++) {
            /* always pull the gaussian from the lower quadrant plane for k = 0
             * plane*/
            /* always pull the whitenoise from the lower quadrant plane for k = 0
             * plane and k == Nmesh / 2 plane*/
            int d1 = 0, d2 = 0;
            int cj = pm->Nmesh[1] - j;
            if(cj >= pm->Nmesh[1]) cj -= pm->Nmesh[1];

            /* d1, d2 points to the conjugate quandrant */
            if( (ci == i && cj < j)
             || (ci < i && cj != j)
             || (ci < i && cj == j)) {
                d1 = 1;
                d2 = 1;
            }

            unsigned int seed_conj, seed_this;
            /* the lower quadrant generator */
            seed_conj = GETSEED(pm, seedtable, i, j, d1, d2);
            gsl_rng_set(lower_rng, seed_conj);

            seed_this = GETSEED(pm, seedtable, i, j, 0, 0);
            gsl_rng_set(this_rng, seed_this);

            for(k = 0; k <= pm->Nmesh[2] / 2; k ++) {
                int use_conj = (d1 != 0 || d2 != 0) && (k == 0 || k == pm->Nmesh[2] / 2);

                double ampl, phase;
                if(use_conj) {
                    /* on k = 0 and Nmesh/2 plane, we use the lower quadrant generator, 
                     * then hermit transform the result if it is nessessary */
                    SAMPLE(this_rng, &ampl, &phase);
                    SAMPLE(lower_rng, &ampl, &phase);
                } else {
                    SAMPLE(lower_rng, &ampl, &phase);
                    SAMPLE(this_rng, &ampl, &phase);
                }

                ptrdiff_t iabs[3] = {i, j, k};
                ptrdiff_t ip = 0;
                for(d = 0; d < 3; d ++) {
                    irel[d] = iabs[d] - pm->ORegion.start[d];
                    ip += pm->ORegion.strides[d] * irel[d];
                }

                if(irel[2] < 0) continue;
                if(irel[2] >= pm->ORegion.size[2]) continue;

                /* we want two numbers that are of std ~ 1/sqrt(2) */
                ampl = sqrt(- log(ampl));

                (delta_k + 2 * ip)[0] = ampl * cos(phase);
                (delta_k + 2 * ip)[1] = ampl * sin(phase);

                if(use_conj) {
                    (delta_k + 2 * ip)[1] *= -1;
                }

                if((pm->Nmesh[0] - iabs[0]) % pm->Nmesh[0] == iabs[0] &&
                   (pm->Nmesh[1] - iabs[1]) % pm->Nmesh[1] == iabs[1] &&
                   (pm->Nmesh[2] - iabs[2]) % pm->Nmesh[2] == iabs[2]) {
                    /* The mode is self conjuguate, thus imaginary mode must be zero */
                    (delta_k + 2 * ip)[1] = 0;
                    (delta_k + 2 * ip)[0] = ampl * cos(phase);
                }

                if(iabs[0] == 0 && iabs[1] == 0 && iabs[2] == 0) {
                    /* the mean is zero */
                    (delta_k + 2 * ip)[0] = 0;
                    (delta_k + 2 * ip)[1] = 0;
                }
            }
        }
        gsl_rng_free(lower_rng);
        gsl_rng_free(this_rng);
    }
    for(i = 0; i < 2; i ++)
    for(j = 0; j < 2; j ++) {
        free(seedtable[i][j]);
    }
/*
    char * fn[1000];
    sprintf(fn, "canvas.dump.f4.%d", pm->ThisTask);
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->ORegion.total * 2, fopen(fn, "w"));
*/
}

static void
pmic_fill_gaussian_fast(PM * pm, FastPMFloat * delta_k, int seed)
{
    ptrdiff_t ind;
    int d;

    gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    /* set uncorrelated seeds */
    gsl_rng_set(random_generator, seed);
    for(d = 0; d < pm->ThisTask * 8; d++) {
        seed = 0x7fffffff * gsl_rng_uniform(random_generator);
    }

    gsl_rng_set(random_generator, seed);

    FastPMFloat * g_x = pm_alloc(pm);

    for(ind = 0; ind < pm->IRegion.total; ind += 2) {
        double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
        double ampl;
        do
            ampl = gsl_rng_uniform(random_generator);
        while(ampl == 0.0);

        /* we need two gaussians of std=1.0 in real space (see footnote 1) */
        ampl = sqrt(-2.0 * log(ampl));
        /* r2c will reduce the variance, so we compensate here. */
        ampl *= sqrt(pm_norm(pm));

        g_x[ind] = ampl * sin(phase);
        g_x[ind + 1] = ampl * cos(phase);
    }
    pm_r2c(pm, g_x, delta_k);
    pm_free(pm, g_x);
}

static void 
pmic_fill_gaussian_slow(PM * pm, FastPMFloat * delta_k, int seed)
{
    ptrdiff_t i[3] = {0};
    int d;
    gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, seed);

    FastPMFloat * g_x = pm_alloc(pm);

    for(i[0] = 0; i[0] < pm->Nmesh[0]; i[0]++)
    for(i[1] = 0; i[1] < pm->Nmesh[1]; i[1]++)
    for(i[2] = 0; i[2] < pm->Nmesh[2]; i[2]++) {
        double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
        double ampl;
        do
            ampl = gsl_rng_uniform(random_generator);
        while(ampl == 0.0);
        ptrdiff_t ii[3];
        ptrdiff_t ind = 0;
        for(d = 0; d < 3; d ++) {
            if(i[d] < pm->IRegion.start[d]) goto next;
            if(i[d] >= pm->IRegion.start[d] + pm->IRegion.size[d]) goto next;
            ii[d] = i[d] - pm->IRegion.start[d];
            ind += ii[d] * pm->IRegion.strides[d];
        }
        /* we need two gaussians of std=1.0 in real space */
        ampl = sqrt(-2.0 * log(ampl));

        /* r2c will reduce the variance, so we compensate here. */
        ampl *= sqrt(pm_norm(pm));

        g_x[ind] = ampl * sin(phase);
        next:
        continue;
    }
    pm_r2c(pm, g_x, delta_k);
    pm_free(pm, g_x);
    gsl_rng_free(random_generator);
}


/* Footnotes */ 

/* 1): 
 * We want delta(k) = delta_real + I delta_imag, where delta_real and
 * delta_imag are Gaussian random variables with variance given by
 * power spectrum, \sigma^2=P(k). We can obtain this equivalently as
 *
 *   delta(k) = A exp(i phase),
 *
 * where the phase is random (i.e. sampled from a uniform distribution)
 * and the amplitude A follows a Rayleigh distribution (see 
 * https://en.wikipedia.org/wiki/Rayleigh_distribution). To sample from 
 * Rayleigh distribution, use inverse transform sampling
 * (see https://en.wikipedia.org/wiki/Inverse_transform_sampling), i.e.
 * start from uniform random variable in [0,1] and then apply inverse of CDF
 * of Rayleigh distribution. From F(A)=CDF(A)=1-e^{-A^2/(2\sigma^2)} we get
 * A = \sigma \sqrt{-2 ln(1-CDF)}. So if x is uniform random number in [0,1], then 
 * A = \sigma \sqrt(-2 ln(x)) follows Rayleigh distribution as desired. 
 * Here we used x instead of 1-x because this does not make a difference for a 
 * uniform random number in [0,1]. In the code below, we start with \sigma=1 and 
 * multiply by sqrt(P(k)) later.
 */
