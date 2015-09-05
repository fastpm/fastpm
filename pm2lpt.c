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
#include "msg.h"

#define PM_2LPT_LOAD_NOISE_K
//#define PM_2LPT_LOAD_DIGRAD
#define PM_2LPT_DUMP

typedef double (*pkfunc)(double k, void * data);

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

    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double k2;
        pm_unravel_o_index(pm, ind, i);
        k2 = index_to_k2(pm, i, k);
        ptrdiff_t i2 = ind * 2;
        /* i k[d] / k2 */
        if(k2 > 0) {
            pm->workspace[i2 + 0] = - pm->canvas[i2 + 1] * (k[dir] / k2);
            pm->workspace[i2 + 1] =   pm->canvas[i2 + 0] * (k[dir] / k2);
        } else {
            pm->workspace[i2 + 0] = 0;
            pm->workspace[i2 + 1] = 0;
        }
#if 0
        if(dir == 0)
            msg_printf(debug, "%ld %ld %ld %ld %g %g %g %g,%g->%g,%g\n", ind, i[0], i[1], i[2], k[0], k[1], k[2],
            pm->canvas[i2][0], pm->canvas[i2][1],
            pm->workspace[i2][0], pm->workspace[i2][1]
            );
#endif
    }
}

static void apply_2lpt_transfer(PM * pm, int dir1, int dir2) {
    ptrdiff_t ind;
    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double k2;
        pm_unravel_o_index(pm, ind, i);
        k2 = index_to_k2(pm, i, k);
        ptrdiff_t i2 = ind * 2;
        if(k2 > 0) {
            pm->workspace[i2 + 0] = pm->canvas[i2 + 0] * (-k[dir1] * k[dir2] / k2);
            pm->workspace[i2 + 1] = pm->canvas[i2 + 1] * (-k[dir1] * k[dir2] / k2);
        } else {
            pm->workspace[i2 + 0] = 0;
            pm->workspace[i2 + 1] = 0;
        }
    }
}

static void pm_2lpt_fill_pdata(PMStore * p, PM * pm) {
    /* fill pdata with a uniform grid, respecting domain given by pm */
    
    ptrdiff_t i[3];
    ptrdiff_t ind;
    int d;
    p->np = 1;
    for(d = 0; d < 3; d++) {
        p->np *= pm->IRegion.size[d];
    }
    if(p->np > p->np_upper) {
        msg_abort(-1, "Need %td particles; %td allocated\n", p->np, p->np_upper);
    }
    ptrdiff_t ptrmax = 0;
    for(ind = 0; ind < pm->allocsize; ind ++){
        ptrdiff_t ptr;
        uint64_t id;
        pm_unravel_i_index(pm, ind, i);
        /* avoid the padded region */
        if(i[2] >= pm->IRegion.size[2]) continue;

        ptr = 0;
        id = 0;
        for(d = 0; d < 3; d ++) {
            ptr = ptr * pm->IRegion.size[d] + i[d];
            id = id * pm->Nmesh[d] + i[d] + pm->IRegion.start[d];
        }
//        msg_aprintf(debug, "Creating particle at offset %td i = %td %td %td\n", ptr, i[0], i[1], i[2]);
        for(d = 0; d < 3; d ++) {
            p->x[ptr][d] = (i[d] + pm->IRegion.start[d]) * (pm->BoxSize[d] / pm->Nmesh[d]);
            p->v[ptr][d] = 0;
            p->dx1[ptr][d] = 0;
            p->dx2[ptr][d] = 0;
            p->acc[d][ptr] = 0;
            p->id[ptr]  = id;
        }
        if(ptr > ptrmax) ptrmax = ptr;
    }
    if(ptrmax + 1 != p->np) {
        msg_abort(-1, "This is an internal error, particle number mismatched with grid. %td != %td\n", ptrmax + 1, p->np);
    }
}

static void pm_2lpt_fill_gaussian(PM * pm, int seed, pkfunc pk, void * pkdata) {
    ptrdiff_t ind;
    int d;

    msg_printf(info, "Filling initial gaussian fluctuations.\n");

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
        /* ensure the fourier space is a normal distribution */
        ampl /= sqrt(pm->Norm);
        ampl *= sqrt((8 * (M_PI * M_PI * M_PI) / pm->Volume));
        pm->canvas[ind] = ampl * sin(phase);
        pm->canvas[ind + 1] = ampl * cos(phase);
    }
#ifdef PM_2LPT_DUMP
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("noise-r.f4", "w"));
#endif
    msg_printf(info, "Transforming to fourier space .\n");
    pm_r2c(pm);
#ifdef PM_2LPT_LOAD_NOISE_K
    fread(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("input-noise-k.f4", "r"));
#endif
#ifdef PM_2LPT_DUMP
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("noise-k.f4", "w"));
#endif
    msg_printf(info, "Inducing correlation.\n");
    msg_printf(debug, "Volume = %g.\n", pm->Volume);
    /* Watch out: PFFT_TRANSPOSED_OUT is y, x, z */
    for(ind = 0; ind < pm->allocsize / 2; ind ++) {
        ptrdiff_t i[3];
        double k[3];
        double knorm = 0;
        pm_unravel_o_index(pm, ind, i);
        double k2 = index_to_k2(pm, i, k);
        knorm = sqrt(k2);
        double f = sqrt(pk(knorm, pkdata));
/*
        msg_printf(debug, "ind = %td, i = %td %td %td, k = %g %g %g, k2 = %g, pk=%g \n",
                ind, i[0], i[1], i[2], k[0], k[1], k[2], k2, pk(knorm, pkdata));
*/
        pm->canvas[2 * ind + 0] *= f;
        pm->canvas[2 * ind + 1] *= f;
    }
#ifdef PM_2LPT_DUMP
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("overdensity-k.f4", "w"));
#endif
}
void pm_2lpt_main(PMStore * p, int Ngrid, double BoxSize, pkfunc pk, int seed, void * pkdata) {
    PM pm;

    PMInit pminit = {
        .Nmesh = Ngrid,
        .BoxSize = BoxSize,
        .NprocX = 0, /* 0 for auto, 1 for slabs */
    };

    pm_pfft_init(&pm, &pminit, &p->iface, MPI_COMM_WORLD);

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

    int DX1[] = {PACK_DX1_X, PACK_DX1_Y, PACK_DX1_Z};
    int DX2[] = {PACK_DX2_X, PACK_DX2_Y, PACK_DX2_Z};
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};
    ptrdiff_t i;
    int d;

    for(d = 0; d < 3; d++) {
        msg_printf(info, "Solving for DX1 axis = %d\n", d);
        apply_za_transfer(&pm, d);

#ifdef PM_2LPT_DUMP
        char * fnames[] = {"dx1-0.f4", "dx1-1.f4", "dx1-2.f4"};
        fwrite(pm.workspace, sizeof(pm.workspace[0]), pm.allocsize, fopen(fnames[d], "w"));
#endif
        pm_c2r(&pm);
        for(i = 0; i < p->np + pgd.nghosts; i ++) {        
            p->dx1[i][d] = pm_readout_one(&pm, p, i);
        }
        pm_reduce_ghosts(&pgd, DX1[d]);
    } 

    real_t * source;
    source = (real_t*) p->iface.malloc(pm.allocsize * sizeof(source[0]));

    memset(source, 0, sizeof(source[0]) * pm.allocsize);

    real_t * field[3];
    for(d = 0; d < 3; d++ )
        field[d] = p->iface.malloc(pm.allocsize * sizeof(field[d][0]));

    for(d = 0; d< 3; d++) {
        msg_printf(info, "Solving for 2LPT axes = %d %d .\n", d, d);
        apply_2lpt_transfer(&pm, d, d);
        pm_c2r(&pm);
        memcpy(field[d], pm.workspace, pm.allocsize * sizeof(field[d][0]));
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
        msg_printf(info, "Solving for 2LPT axes = %d %d .\n", d1, d2);
        apply_2lpt_transfer(&pm, d1, d2);
        pm_c2r(&pm);
        for(i = 0; i < pm.allocsize; i ++) {
            source[i] -= pm.workspace[i] * pm.workspace[i];
        }
    } 
    memcpy(pm.canvas, source, pm.allocsize * sizeof(pm.canvas[0]));

#ifdef PM_2LPT_DUMP
    fwrite(pm.canvas, sizeof(pm.canvas[0]), pm.allocsize, fopen("digrad.f4", "w"));
#endif
#ifdef PM_2LPT_LOAD_DIGRAD
    fread(pm.canvas, sizeof(pm.canvas[0]), pm.allocsize, fopen("input-digrad.f4", "r"));
#endif
    pm_r2c(&pm);

    for(d = 2; d >=0; d-- )
        p->iface.free(field[d]);
    p->iface.free(source);

    for(d = 0; d < 3; d++) {
        msg_printf(info, "Solving for DX2 axis = %d .\n", d);
        /* 
         * We absorb some the negative factor in za transfer to below;
         *
         * */
        apply_za_transfer(&pm, d);
#ifdef PM_2LPT_DUMP
        char * fnames[] = {"dx2-0.f4", "dx2-1.f4", "dx2-2.f4"};
        fwrite(pm.workspace, sizeof(pm.workspace[0]), pm.allocsize, fopen(fnames[d], "w"));
#endif
        pm_c2r(&pm);
        for(i = 0; i < p->np + pgd.nghosts; i ++) {        
            /* this ensures x = x0 + dx1(t) + 3/ 7 dx2(t) */
            p->dx2[i][d] = 3.0 / 7.0 * pm_readout_one(&pm, p, i) / pm.Norm ;
        }
        pm_reduce_ghosts(&pgd, DX2[d]);
    }
    double dx1disp[3] = {0};
    double dx2disp[3] = {0};

    for(i = 0; i < p->np; i ++) {
        for(d =0; d < 3; d++) {
            dx1disp[d] += p->dx1[i][d] * p->dx1[i][d];
            dx2disp[d] += p->dx2[i][d] * p->dx2[i][d];
        } 
    }
    uint64_t Ntot = p->np;
    MPI_Allreduce(MPI_IN_PLACE, dx1disp, 3, MPI_DOUBLE, MPI_SUM, pm.Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, dx2disp, 3, MPI_DOUBLE, MPI_SUM, pm.Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, &Ntot,   1, MPI_LONG,  MPI_SUM, pm.Comm2D);
    for(d =0; d < 3; d++) {
        dx1disp[d] /= Ntot;
        dx1disp[d] = sqrt(dx1disp[d]);
        dx2disp[d] /= Ntot;
        dx2disp[d] = sqrt(dx2disp[d]);
    }
    msg_printf(info, "dx1 disp : %g %g %g %g\n", 
            dx1disp[0], dx1disp[1], dx1disp[2],
            (dx1disp[0] + dx1disp[1] + dx1disp[2]) / 3.0);
    msg_printf(info, "dx2 disp : %g %g %g %g\n", 
            dx2disp[0], dx2disp[1], dx2disp[2],
            (dx2disp[0] + dx2disp[1] + dx2disp[2]) / 3.0);

#ifdef PM_2LPT_DUMP
    fwrite(p->dx1, sizeof(p->dx1[0]), p->np, fopen("dx1.f4x3", "w"));
    fwrite(p->dx2, sizeof(p->dx2[0]), p->np, fopen("dx2.f4x3", "w"));
#endif
}

