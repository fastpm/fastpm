#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>

#include "libfastpm.h"
#include "pmpfft.h"
#include "vpm.h"
#include "msg.h"

VPM *
vpm_find(VPM * vpm, double a) 
{
    /* find the PM object for force calculation at time a*/
    int i;
    for (i = 0; !vpm[i].end; i ++) {
        if(vpm[i].a_start > a) break;
    }
    return &vpm[i-1];
}

VPM * 
vpm_create (int size, int * pm_nc_factors, double * change_pm, 
        PMInit * baseinit, PMIFace * iface, MPI_Comm comm) 
{
    /* plan for the variable PMs; keep in mind we do variable
     * mesh size (PM resolution). We plan them at the begining of the run
     * in theory we can use some really weird PM resolution as function
     * of time, but the parameter files / API need to support this.
     * */
    VPM * vpm = malloc(sizeof(VPM) * (size + 1));
    int i;
    for (i = 0; i < size; i ++) {
        vpm[i].end = 0; /* not the last item */

        vpm[i].pm_nc_factor = pm_nc_factors[i];
        vpm[i].a_start = change_pm[i];

        PMInit pminit = *baseinit;
        pminit.Nmesh = baseinit->Nmesh * vpm[i].pm_nc_factor;
        msg_printf(debug, "Nmesh = %td at a %5.4g \n", pminit.Nmesh, change_pm[i]);
        pm_init(&vpm[i].pm, &pminit, iface, comm);
    }
    /* the end of the list */
    vpm[i].end = 1;
    return vpm;
}
void vpm_free (VPM * vpm) {
    int i;
    while(!vpm[i].end) {
        pm_destroy(&vpm[i].pm);

        i++;
    }
    free(vpm);
}
#if 0
double
vpm_estimate_alloc_factor(double failure_rate) 
{
    double factor = 1.0;
    int i;
    for(i = 0; i < n_vpm; i ++) {
        PM * pm = &_vpm[i].pm;
        double Volume = pm->Volume / pm->NTask;

        /* Correct for the surface area */
        Volume += 2 * pm->Volume / pm->Nmesh[0] * pm->Nproc[0];
        Volume += 2 * pm->Volume / pm->Nmesh[1] * pm->Nproc[1];

        double factor1 = estimate_alloc_factor(Volume, pm->NTask, failure_rate);
        if(factor1 > factor) factor = factor1;
    }
    return factor;
}
#endif
