#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmstore.h"
#include "pmpfft.h"
#include "pmghosts.h"


FastPMFloat * pm_alloc(PM * pm) 
{
    return pm->iface.malloc(sizeof(FastPMFloat) * pm->allocsize);
}

void 
pm_free(PM * pm, FastPMFloat * data) 
{
    pm->iface.free(data);
}

void 
pm_assign(PM * pm, FastPMFloat * from, FastPMFloat * to) 
{
    memcpy(to, from, sizeof(from[0]) * pm->allocsize);
}

size_t 
pm_size(PM * pm)
{
    return pm->allocsize;
}

ptrdiff_t * pm_nmesh(PM * pm) {
    return pm->Nmesh;
}
double * pm_boxsize(PM * pm) {
    return pm->BoxSize;
}

PMRegion * pm_i_region(PM * pm) {
    return &pm->IRegion;
}

PMRegion * pm_o_region(PM * pm) {
    return &pm->ORegion;
}

void
fastpm_utils_paint(PM * pm, PMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k) 
{

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    /* since for 2lpt we have on average 1 particle per cell, use 1.0 here.
     * otherwise increase this to (Nmesh / Ngrid) **3 */
    FastPMFloat * canvas = pm_alloc(pm);
    pm_paint(pm, canvas, p, p->np + pgd->nghosts, 1.0);
    pm_assign(pm, canvas, delta_x);

    if(delta_k) {
        pm_r2c(pm, canvas, delta_k);
        ptrdiff_t i = 0;
        for(i = 0; i < pm->allocsize; i ++) {
            delta_k[i] /= pm->Norm;
        }
    }
    pm_free(pm, canvas);
    pm_ghosts_free(pgd);
}

void 
fastpm_utils_dump(PM * pm , char * filename, FastPMFloat *data) 
{
    char fn1[1024];
    char fn2[1024];
    if(pm->NTask > 1) {
        sprintf(fn1, "%s.%03d.geometry", filename, pm->ThisTask);
        sprintf(fn2, "%s.%03d", filename, pm->ThisTask);
    } else {
        sprintf(fn1, "%s.geometry", filename);
        sprintf(fn2, "%s", filename);
    }
    FILE * fp;
    fp = fopen(fn2, "w");
    fwrite(data, sizeof(FastPMFloat), pm->allocsize, fp);
    fclose(fp);
    fp = fopen(fn1, "w");
    fprintf(fp, "%td %td %td\n", 
                    pm->IRegion.strides[0],
                    pm->IRegion.strides[1],
                    pm->IRegion.strides[2]);
    fprintf(fp, "%td %td %td\n", 
                    pm->IRegion.size[0],
                    pm->IRegion.size[1],
                    pm->IRegion.size[2]);
    fclose(fp);
}

static double 
tk_eh(double k, struct fastpm_powerspec_eh_params * params)		/* from Martin White */
{
    double q, theta, ommh2, a, s, gamma, L0, C0;
    double tmp;
    double omegam, ombh2, hubble;

    /* other input parameters */
    hubble = params->hubble_param;
    omegam = params->omegam;
    ombh2 = params->omegab * hubble * hubble;

    theta = 2.728 / 2.7;
    ommh2 = omegam * hubble * hubble;
    s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
    a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
        + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
    gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
    gamma *= omegam * hubble;
    q = k * theta * theta / gamma;
    L0 = log(2. * exp(1.) + 1.8 * q);
    C0 = 14.2 + 731. / (1. + 62.5 * q);
    tmp = L0 / (L0 + C0 * q * q);
    return (tmp);
}

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param)	/* Eisenstein & Hu */
{
    return param->Norm * k * pow(tk_eh(k, param), 2);
}

