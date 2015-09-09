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
#include "pm2lpt.h"
#include "msg.h"

//#define PM_2LPT_LOAD_NOISE_K
//#define PM_2LPT_LOAD_DIGRAD
//#define PM_2LPT_DUMP

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

static void 
apply_za_transfer(PM * pm, int dir) 
{
    /* apply za transfer function i k / k2 from canvas -> workspace */
    ptrdiff_t ind;
    ptrdiff_t i[3] = {0};
    for(ind = 0; ind < pm->ORegion.total * 2; ind += 2, pm_inc_o_index(pm, i)) {
        double k[3];
        double k2;
        k2 = index_to_k2(pm, i, k);
        /* i k[d] / k2 */
        if(k2 > 0) {
            pm->workspace[ind + 0] = - pm->canvas[ind + 1] * (k[dir] / k2);
            pm->workspace[ind + 1] =   pm->canvas[ind + 0] * (k[dir] / k2);
        } else {
            pm->workspace[ind + 0] = 0;
            pm->workspace[ind + 1] = 0;
        }
#if 0
        if(dir == 0)
            msg_printf(debug, "%ld %ld %ld %ld %g %g %g %g,%g->%g,%g\n", ind, i[0], i[1], i[2], k[0], k[1], k[2],
            pm->canvas[ind][0], pm->canvas[ind][1],
            pm->workspace[ind][0], pm->workspace[ind][1]
            );
#endif
    }
}

static void 
apply_2lpt_transfer(PM * pm, int dir1, int dir2) 
{
    /* apply 2lpt transfer function - k k / k2 from canvas -> workspace */
    ptrdiff_t ind;
    ptrdiff_t i[3] = {0};
    for(ind = 0; ind < pm->ORegion.total * 2; ind += 2, pm_inc_o_index(pm, i)) {
        double k[3];
        double k2;
        k2 = index_to_k2(pm, i, k);
        if(k2 > 0) {
            pm->workspace[ind + 0] = pm->canvas[ind + 0] * (-k[dir1] * k[dir2] / k2);
            pm->workspace[ind + 1] = pm->canvas[ind + 1] * (-k[dir1] * k[dir2] / k2);
        } else {
            pm->workspace[ind + 0] = 0;
            pm->workspace[ind + 1] = 0;
        }
    }
}

static void 
pm_2lpt_fill_pdata(PMStore * p, PM * pm) 
{
    /* fill pdata with a uniform grid, respecting domain given by pm */
    
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
    ptrdiff_t i[3] = {0};
    ptrdiff_t ptr = 0;
    for(ind = 0; ind < pm->IRegion.total; ind ++, pm_inc_i_index(pm, i)){
        uint64_t id;
        /* avoid the padded region */
        if(i[2] >= pm->IRegion.size[2]) continue;

        id = 0;
        for(d = 0; d < 3; d ++) {
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
        ptr ++;
    }
    if(ptrmax + 1 != p->np) {
        msg_abort(-1, "This is an internal error, particle number mismatched with grid. %td != %td, allocsize=%td, shape=(%td %td %td)\n", 
            ptrmax + 1, p->np, pm->allocsize, 
            pm->IRegion.size[0],
            pm->IRegion.size[1],
            pm->IRegion.size[2]

            );
    }
}
inline void 
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
inline unsigned int 
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
pm_2lpt_fill_gaussian_gadget(PM * pm, int seed, pkfunc pk, void * pkdata) 
{
    /* Fill gaussian with gadget scheme */
    ptrdiff_t ind;
    int d;
    int i, j, k;

    memset(pm->canvas, 0, sizeof(pm->canvas[0]) * pm->allocsize);

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
                pm->canvas[2 * ip + 0] = delta * cos(phase);
                pm->canvas[2 * ip + 1] = delta * sin(phase);

                if(hermitian && k == 0) {
                    pm->canvas[2 * ip + 1] *= -1;
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
pm_2lpt_fill_gaussian(PM * pm, int seed, pkfunc pk, void * pkdata) 
{
    ptrdiff_t ind;
    int d;

    msg_printf(info, "Filling initial white noise field in x-space.\n");
    gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, seed * pm->NTask + pm->ThisTask);
    ptrdiff_t i[3] = {0};

    /* first fill in x-space, then transform to ensure conjugates;
     * this means no zoomins! */
    for(ind = 0; ind < pm->IRegion.total; ind += 2) {
        double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
        double ampl;
        do
            ampl = gsl_rng_uniform(random_generator);
        while(ampl == 0.0);

        ampl = sqrt(-log(ampl));
        /* ensure the fourier space is a normal distribution */
        ampl /= sqrt(pm->Norm);
        /* 2pi / k -- this matches the dimention of sqrt(p) but I always 
         * forget where it is from. */
        ampl *= sqrt((8 * (M_PI * M_PI * M_PI) / pm->Volume));
        pm->workspace[ind] = ampl * sin(phase);
        pm->workspace[ind + 1] = ampl * cos(phase);
    }
    
    gsl_rng_set(random_generator, seed);
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
        ampl /= sqrt(pm->Norm);

        ampl *= sqrt((8 * (M_PI * M_PI * M_PI) / pm->Volume));
        pm->workspace[ind] = ampl * sin(phase);
        next:
        continue;
    }

#ifdef PM_2LPT_DUMP
    fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen("noise-r.f4", "w"));
#endif
    msg_printf(info, "Transforming to fourier space .\n");
    pm_r2c(pm);

#ifdef PM_2LPT_LOAD_NOISE_K
    fread(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("input-noise-k.f4", "r"));
#endif
#ifdef PM_2LPT_DUMP
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("noise-k.f4", "w"));
#endif
    msg_printf(info, "Inducing correlation to the white noise.\n");
    msg_printf(debug, "Volume = %g.\n", pm->Volume);

    i[0] = i[1] = i[2] = 0;
    for(ind = 0; ind < pm->ORegion.total *2 ; ind += 2, pm_inc_o_index(pm, i)) {
        double k[3];
        double knorm = 0;
        double k2 = index_to_k2(pm, i, k);
        knorm = sqrt(k2);
        double f = sqrt(pk(knorm, pkdata));

//        msg_aprintf(debug, "ind = %td, i = %td %td %td, k = %g %g %g, k2 = %g, pk=%g \n",
//                ind, i[0], i[1], i[2], k[0], k[1], k[2], k2, pk(knorm, pkdata));

        pm->canvas[ind + 0] *= f;
        pm->canvas[ind + 1] *= f;
    }

#ifdef PM_2LPT_DUMP
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("overdensity-k.f4", "w"));
#endif
}

void 
pm_2lpt_main(PMStore * p, int Ngrid, double BoxSize, pkfunc pk, 
        int seed, void * pkdata, MPI_Comm comm) 
{
    PM pm;

    PMInit pminit = {
        .Nmesh = Ngrid,
        .BoxSize = BoxSize,
        .NprocX = 0, /* 0 for auto, 1 for slabs */
        .transposed = 0,
    };

    pm_pfft_init(&pm, &pminit, &p->iface, comm);

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
    
    pm_2lpt_fill_gaussian_gadget(&pm, seed, pk, pkdata);
//    pm_2lpt_fill_gaussian(&pm, seed, pk, pkdata);

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
        for(i = 0; i < pm.IRegion.total; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }    
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
        msg_printf(info, "Solving for 2LPT axes = %d %d .\n", d1, d2);
        apply_2lpt_transfer(&pm, d1, d2);
        pm_c2r(&pm);
        for(i = 0; i < pm.IRegion.total; i ++) {
            source[i] -= pm.workspace[i] * pm.workspace[i];
        }
    } 
    memcpy(pm.workspace, source, pm.allocsize * sizeof(pm.canvas[0]));

#ifdef PM_2LPT_DUMP
    fwrite(pm.workspace, sizeof(pm.workspace[0]), pm.allocsize, fopen("digrad.f4", "w"));
#endif
#ifdef PM_2LPT_LOAD_DIGRAD
    fread(pm.workspace, sizeof(pm.workspace[0]), pm.allocsize, fopen("input-digrad.f4", "r"));
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
            p->dx2[i][d] = pm_readout_one(&pm, p, i) / pm.Norm ;
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

