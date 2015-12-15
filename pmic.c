#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <gsl/gsl_rng.h>

#include "pmpfft.h"
#include "pmic.h"
#include "msg.h"

static double 
index_to_k2(PM * pm, ptrdiff_t i[3], double k[3]) 
{
    /* convert (rel) integer to k values. */
    int d;
    double k2 = 0;
    for(d = 0; d < 3 ; d++) {
        k[d] = pm->MeshtoK[d][i[d] + pm->ORegion.start[d]];
        k2 += k[d] * k[d];
    } 
    return k2;
}


static inline void 
SETSEED(PM * pm, unsigned int * table[2][2], int i, int j, gsl_rng * rng) 
{ 
    unsigned int seed = 0x7fffffff * gsl_rng_uniform(rng); 

    int ii[2] = {i, pm->Nmesh[0] - i};
    int jj[2] = {j, pm->Nmesh[1] - j};
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

void 
pmic_fill_gaussian_gadget(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata) 
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

    double fac = sqrt(pow(2 * M_PI, 3) / pm->Volume);

    ptrdiff_t irel[3];
    for(i = pm->ORegion.start[0]; 
        i < pm->ORegion.start[0] + pm->ORegion.size[0]; 
        i ++) {

        gsl_rng * random_generator0 = gsl_rng_alloc(gsl_rng_ranlxd1);
        gsl_rng * random_generator1 = gsl_rng_alloc(gsl_rng_ranlxd1);

        for(j = pm->ORegion.start[1]; 
            j < pm->ORegion.start[1] + pm->ORegion.size[1]; 
            j ++) {
            /* always pull the gaussian from the lower quadrant plane for k = 0
             * plane*/
            int hermitian = 0;
            int d1 = 0, d2 = 0;

            if(i == 0) {
                if(j > pm->Nmesh[1] / 2) {
                    hermitian = 1; 
                    d2 = 1;
                } 
            } else {
                if(i > pm->Nmesh[0] / 2) {
                    hermitian = 1;
                    d1 = 1;
                    d2 = j != 0;
                }  else {
                    /* no transpose */
                    d1 = d2 = 0;
                }
            } 

            unsigned int seed;
            seed = GETSEED(pm, seedtable, i, j, d1, d2);
            gsl_rng_set(random_generator0, seed);

            seed = GETSEED(pm, seedtable, i, j, 0, 0);
            gsl_rng_set(random_generator1, seed);

            /* this black magic matches two generators. */
            double skip = gsl_rng_uniform(random_generator1);
            do skip = gsl_rng_uniform(random_generator1);
            while(skip == 0);

            for(k = 0; k <= pm->Nmesh[2] / 2; k ++) {
                gsl_rng * random_generator = k?random_generator1:random_generator0;
                /* on k = 0 plane, we use the lower quadrant generator, 
                 * then hermit transform the result if it is nessessary */
                double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
                double ampl = 0;
                do ampl = gsl_rng_uniform(random_generator); while(ampl == 0);

                ptrdiff_t iabs[3] = {i, j, k};
                ptrdiff_t ip = 0;
                for(d = 0; d < 3; d ++) {
                    irel[d] = iabs[d] - pm->ORegion.start[d];
                    ip += pm->ORegion.strides[d] * irel[d];
                }

                if(irel[2] < 0) continue;
                if(irel[2] >= pm->ORegion.size[2]) continue;

                double tmp[3];
                double k2 = index_to_k2(pm, irel, tmp);
                double kmag = sqrt(k2);
                double p_of_k = - log(ampl);

                for(d = 0; d < 3; d ++) {
                    if(iabs[d] == pm->Nmesh[d] / 2) p_of_k = 0;
                } 
                if(iabs[0] == 0 && iabs[1] == 0 && iabs[2] == 0) {
                    p_of_k = 0;
                }
            
                p_of_k *= pk(kmag, pkdata);

                double delta = fac * sqrt(p_of_k);
                delta_k[2 * ip + 0] = delta * cos(phase);
                delta_k[2 * ip + 1] = delta * sin(phase);

                if(hermitian && k == 0) {
                    delta_k[2 * ip + 1] *= -1;
                }
            }
        }
        gsl_rng_free(random_generator0);
        gsl_rng_free(random_generator1);
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
pmic_induce_correlation(PM * pm, float_t * g_x, float_t * delta_k, pkfunc pk, void * pkdata) {

    msg_printf(info, "Transforming to fourier space .\n");
    pm_r2c(pm, g_x, delta_k);

    msg_printf(info, "Inducing correlation to the white noise.\n");
    msg_printf(debug, "Volume = %g.\n", pm->Volume);

    ptrdiff_t i[3] = {0};
    ptrdiff_t ind;
    for(ind = 0; ind < pm->ORegion.total *2 ; ind += 2, pm_inc_o_index(pm, i)) {
        double k[3];
        double knorm = 0;
        double k2 = index_to_k2(pm, i, k);
        knorm = sqrt(k2);
        double f = sqrt(pk(knorm, pkdata));

        /* ensure the fourier space is a normal distribution */
        f /= sqrt(pm->Norm);
        /* 2pi / k -- this matches the dimention of sqrt(p) but I always 
         * forget where it is from. */
        f *= sqrt((8 * (M_PI * M_PI * M_PI) / pm->Volume));
//
//        msg_aprintf(debug, "ind = %td, i = %td %td %td, k = %g %g %g, k2 = %g, pk=%g \n",
//                ind, i[0], i[1], i[2], k[0], k[1], k[2], k2, pk(knorm, pkdata));

        delta_k[ind + 0] *= f;
        delta_k[ind + 1] *= f;
    }
}


void 
pmic_fill_gaussian_fast(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata)
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

    float_t * g_x = pm_alloc(pm);

    for(ind = 0; ind < pm->IRegion.total; ind += 2) {
        double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
        double ampl;
        do
            ampl = gsl_rng_uniform(random_generator);
        while(ampl == 0.0);

        ampl = sqrt(-log(ampl));
        g_x[ind] = ampl * sin(phase);
        g_x[ind + 1] = ampl * cos(phase);
    }
    pmic_induce_correlation(pm, g_x, delta_k, pk, pkdata);
    pm_free(pm, g_x);
}

void 
pmic_fill_gaussian_slow(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata) 
{    
    ptrdiff_t i[3] = {0};
    int d;
    gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, seed);

    float_t * g_x = pm_alloc(pm);

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
        ampl = sqrt(-log(ampl));
        g_x[ind] = ampl * sin(phase);
        next:
        continue;
    }
    pmic_induce_correlation(pm, g_x, delta_k, pk, pkdata);
    pm_free(pm, g_x);
    gsl_rng_free(random_generator);
}
            
void 
pmic_read_gaussian(PM * pm, float_t * delta_k, char * filename, pkfunc pk, void * pkdata)
{
    ptrdiff_t ind;
    int d;
    ptrdiff_t i[3] = {0, 0, 0};

    msg_printf(info, "Reading grafic white noise file from '%s'.\n", filename);

    FILE * fp = fopen(filename, "r");

    struct {
        int32_t bs1;
        int32_t n[3];
        int32_t seed;
        int32_t bs2;
    } header;

    fread(&header, sizeof(header), 1, fp);

    ind = 0;

    if(header.bs1 != 16) 
        msg_abort(-1, "file not in BigMD noise format.\n");

    msg_printf(info, "GrafIC noise is Fortran ordering. FastPM is in C ordering."
        "The simulation will be transformed x->z y->y z->x.\n");

    for(d = 0; d < 3; d ++) {
        /* BigMD is in */
        int permute[3] = {2, 1, 0};
        if(header.n[d] != pm->Nmesh[permute[d]]) {
            msg_abort(-1, "file is in %d, but simulation is in %d.\n",
                header.n[d], pm->Nmesh[permute[d]]);
        }
    }

    float * buf = malloc(sizeof(float) * header.n[0] * pm->IRegion.size[1]);

    float_t * g_x = pm_alloc(pm);

    for(i[0] = 0; i[0] < pm->IRegion.size[0]; i[0] ++) {
        ptrdiff_t i_abs[3];

        i_abs[0] = i[0] + pm->IRegion.start[0];

        ptrdiff_t offset = (/*header*/16 + 8) + 
                    i_abs[0] * (header.n[0] * header.n[1] * 4 + 8); 
        int32_t bs;

        fseek(fp, offset, SEEK_SET);

        fread(&bs, 4, 1, fp);

        if(bs != 4 * header.n[0] * header.n[1])
            msg_abort(-1, "file size is wrong\n");

        fseek(fp, pm->IRegion.start[1] * header.n[1] * 4, SEEK_CUR);
        fread(buf, 4, header.n[0] * pm->IRegion.size[1], fp);
        int p = 0;
        for(i[1] = 0; i[1] < pm->IRegion.size[1]; i[1] ++) {
            /* note that size[2] is Nmesh[2] == n[0] and start[2] is 0 */
            for(i[2] = 0; i[2] < pm->IRegion.size[2]; i[2] ++) {
                ind = 0;
                for(d = 0; d < 3; d++) {
                    ind += i[d] * pm->IRegion.strides[d];
                }
                g_x[ind] = buf[p];
                p ++;
            }
        }
    }
    free(buf);
    MPI_Barrier(pm->Comm2D);    

    pmic_induce_correlation(pm, g_x, delta_k, pk, pkdata);
    pm_free(pm, g_x);
}

