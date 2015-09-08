#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include "pmpfft.h"

/* paint and read out */
static void pm_paint_pos(PM * pm, double pos[3], double weight) {
    /* to pm->workspace */
    int n;
    double gpos[3];
    int ipos[3];
    double k[3][2];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->Nmesh[d] / pm->BoxSize[d];
        ipos[d] = floor(gpos[d]);
        k[d][0] = 1 + ipos[d] - gpos[d];
        k[d][1] = gpos[d] - ipos[d];
    }
    for(n = 0; n < 8; n ++) {
        double kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n>>d) & 1;

            kernel *= k[d][rel];

            int targetpos = ipos[d] + rel - pm->IRegion.start[d];

            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) 
                goto outside;
            ind += pm->IRegion.strides[d] * targetpos;
        }
        pm->workspace[ind] += weight * kernel;

    outside:
        continue;
    }
    return;
}

static double pm_readout_pos(PM * pm, double pos[3]) {
    /* from pm->workspace */
    double value = 0;
    int n;
    double gpos[3];
    int ipos[3];
    double k[3][2];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->Nmesh[d] / pm->BoxSize[d];
        ipos[d] = floor(gpos[d]);
        k[d][0] = 1 + ipos[d] - gpos[d];
        k[d][1] = gpos[d] - ipos[d];
    }
    for(n = 0; n < 8; n ++) {
        double  kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n>>d) & 1;

            kernel *= k[d][rel];

            int targetpos = ipos[d] + rel - pm->IRegion.start[d];

            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) {
                goto outside;
            } 
            ind += pm->IRegion.strides[d] * targetpos;
        }
        value += kernel * pm->workspace[ind];
outside:
        continue;
    }
    return value;
}

void pm_paint(PM * pm, void * pdata, ptrdiff_t size) {
    ptrdiff_t i;
    memset(pm->workspace, 0, sizeof(pm->workspace[0]) * pm->allocsize);
    for (i = 0; i < size; i ++) {
        double pos[3];
        pm->iface.get_position(pdata, i, pos);

/*
        if(pm_pos_to_rank(pm, pos) != pm->ThisTask) {
            raise(SIGTRAP);
        }
*/
        pm_paint_pos(pm, pos, 1.0);
    }
}

double pm_readout_one(PM * pm, void * pdata, ptrdiff_t i) {
    double pos[3];
    pm->iface.get_position(pdata, i, pos);
    return pm_readout_pos(pm, pos);
}
