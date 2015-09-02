#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>

#include "pmpfft.h"

/* paint and read out */
static void pm_paint_pos(PM * pm, double pos[3], double weight) {
    /* to pm->canvas */
    int n;
    for(n = 0; n < 8; n ++) {
        double  kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            double gpos = pos[d] * pm->Nmesh[d] / pm->BoxSize[d];
            int intpos = floor(gpos);
            double diff = gpos - intpos;
            int rel = (n>>d) & 1;
            int targetpos = intpos + rel;
            if(rel != 0) {
                kernel *= diff;
            } else {
                kernel *= (1.0 - diff);
            }
            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) {
                goto outside;
            } 
            ind = ind * pm->IRegion.size[d] + targetpos;
        }
        pm->canvas[ind] += weight * kernel;
    }
outside:
    return;
}

static double pm_readout_pos(PM * pm, double pos[3]) {
    /* from pm->workspace */
    double value = 0;
    int n;
    for(n = 0; n < 8; n ++) {
        double  kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            double gpos = pos[d] * pm->Nmesh[d] / pm->BoxSize[d];
            int intpos = floor(gpos);
            double diff = gpos - intpos;
            int rel = (n>>d) & 1;
            int targetpos = intpos + rel;
            if(rel != 0) {
                kernel *= diff;
            } else {
                kernel *= (1.0 - diff);
            }
            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) {
                goto outside;
            } 
            ind = ind * pm->IRegion.size[d] + targetpos;
        }
        value += kernel * pm->workspace[ind];
    }
    return value;
outside:
    return 0;
}

void pm_paint(PM * pm, void * pdata, ptrdiff_t size) {
    ptrdiff_t i;
    memset(pm->canvas, 0, sizeof(pm->canvas[0]) * pm->allocsize);
    for (i = 0; i < size; i ++) {
        double pos[3];
        pm->iface.get_position(pdata, i, pos);
        pm_paint_pos(pm, pos, 1.0);
    }
}

double pm_readout_one(PM * pm, void * pdata, ptrdiff_t i) {
    double pos[3];
    pm->iface.get_position(pdata, i, pos);
    return pm_readout_pos(pm, pos);
}
