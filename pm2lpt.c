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

typedef double (*pkfunc)(double k, void * data);

static void unravel_o_index(PM * pm, ptrdiff_t ind, ptrdiff_t i[3]) {
    ptrdiff_t tmp = ind;
    i[1] = tmp / pm->ORegion.strides[1];
    tmp %= pm->ORegion.strides[1];
    i[0] = tmp / pm->ORegion.strides[0];
    tmp %= pm->ORegion.strides[0];
    i[2] = tmp;
}
static double index_to_k2(PM * pm, ptrdiff_t i[3], double k[3]) {
    int d;
    double k2 = 0;
    for(d = 0; d < 3 ; d++) {
        k[d] = pm->MeshtoK[d][i[d] + pm->ORegion.start[d]];
        k2 += k[d] * k[d];
    } 
    return k2;
}

static void apply_za_transfer(PM * pm, int dir) {
    ptrdiff_t ind;
    pfft_complex * c1 = (pfft_complex*) pm->canvas;
    pfft_complex * c2 = (pfft_complex*) pm->workspace;

    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double k2;
        unravel_o_index(pm, ind, i);
        k2 = index_to_k2(pm, i, k);

        /* i k[d] / k2 */
        c2[ind][0] = - c1[ind][1] * (k[dir] / k2);
        c2[ind][1] =   c1[ind][0] * (k[dir] / k2);
    }
}

static void apply_2lpt_transfer(PM * pm, int dir1, int dir2) {
    ptrdiff_t ind;
    pfft_complex * c1 = (pfft_complex*) pm->canvas;
    pfft_complex * c2 = (pfft_complex*) pm->workspace;

    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double k2;
        unravel_o_index(pm, ind, i);
        k2 = index_to_k2(pm, i, k);

        c2[ind][0] = c1[ind][0] * (-k[dir1] * k[dir2] / k2);
        c2[ind][1] = c1[ind][1] * (-k[dir1] * k[dir2] / k2);
    }
}

static void pm_2lpt_fill_pdata(PMStore * p, PM * pm) {
    /* fill pdata with a uniform grid, respecting domain given by pm */
    pm_store_init(p, 10);
    p->np = 0;
}

static void pm_2lpt_fill_gaussian(PM * pm, int seed, pkfunc pk, void * pkdata) {
    ptrdiff_t ind;
    int d;
    gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, seed * pm->NTask + pm->ThisTask);

    /* first fill in x-space, then transform to ensure conjugates;
     * this means no zoomins! */
    for(ind = 0; ind < pm->allocsize; ind += 2) {
        double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
        double ampl;
        do
            ampl = gsl_rng_uniform(random_generator);
        while(ampl == 0.0);

        ampl = -log(ampl);
        pm->canvas[ind] = ampl * sin(phase);
        pm->canvas[ind + 1] = ampl * cos(phase);
    }

    pm_r2c(pm);

    /* Watch out: PFFT_TRANSPOSED_OUT is y, x, z */
    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double knorm = 0;
        unravel_o_index(pm, ind, i);
        for(d = 0; d < 3; d ++) {
            k[d] = pm->MeshtoK[d][i[d]];
            knorm += k[d] * k[d];
        } 
        knorm = sqrt(knorm);
        pm->canvas[ind] *= pk(knorm, pkdata);
    }
}
void pm_2lpt_main(PMStore * p, int Ngrid, double BoxSize, pkfunc pk, int seed, void * pkdata) {
    PM pm;

    PMInit pminit = {
        .Nmesh = Ngrid,
        .BoxSize = BoxSize,
        .NprocX = 0, /* 0 for auto, 1 for slabs */
    };

    pm_pfft_init(&pm, &pminit, &p->iface, MPI_COMM_WORLD);

    int DX1[] = {PACK_DX1_X, PACK_DX1_Y, PACK_DX1_Z};
    int DX2[] = {PACK_DX2_X, PACK_DX2_Y, PACK_DX2_Z};
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};
    ptrdiff_t i;
    int d;
    pm_2lpt_fill_pdata(p, &pm);
    
    PMGhostData pgd = {
        .pm = &pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
    };

    pm_append_ghosts(&pgd);

    pm_start(&pm);
    
    pm_2lpt_fill_gaussian(&pm, seed, pk, pkdata);

    for(d = 0; d < 3; d++) {
        apply_za_transfer(&pm, d);
        pm_c2r(&pm);
        for(i = 0; i < p->np + pgd.nghosts; i ++) {        
            pm_readout_one(&pm, p, i);
        }
        pm_reduce_ghosts(&pgd, DX1[d]);
    } 

    double * source;
    source = (double*) p->iface.malloc(pm.allocsize * sizeof(double));

    memset(source, 0, sizeof(double) * pm.allocsize);

    double * field[3];
    for(d = 0; d < 3; d++ )
        field[d] = p->iface.malloc(pm.allocsize * sizeof(double));

    for(d = 0; d< 3; d++) {
        apply_2lpt_transfer(&pm, d, d);
        pm_c2r(&pm);
        memcpy(field[d], pm.workspace, pm.allocsize * sizeof(double));
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
        for(i = 0; i < pm.allocsize; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }    
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
        apply_2lpt_transfer(&pm, d1, d2);
        pm_c2r(&pm);
        for(i = 0; i < pm.allocsize; i ++) {
            source[i] -= pm.workspace[i] * pm.workspace[i];
        }
    } 
    memcpy(pm.canvas, source, pm.allocsize * sizeof(double));

    for(d = 2; d >=0; d-- )
        p->iface.free(field[d]);
    p->iface.free(source);

    for(d = 0; d < 3; d++) {
        apply_za_transfer(&pm, d);
        pm_c2r(&pm);
        for(i = 0; i < p->np + pgd.nghosts; i ++) {        
            pm_readout_one(&pm, p, i);
        }
        pm_reduce_ghosts(&pgd, DX2[d]);
    }
}

