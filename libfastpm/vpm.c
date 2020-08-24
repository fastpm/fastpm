#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"
#include "vpm.h"

VPM *
vpm_find(VPM * vpm, double a) 
{
    /* find the PM object for force calculation at time a*/
    int i;
    for (i = 0; !vpm[i].end; i ++) {
        if(vpm[i].a_start > a) break;
    }
    /* Start with the first pm. */
    if (i == 0) i = 1;
    return &vpm[i-1];
}

VPM * 
vpm_create (VPMInit * vpminit, PMInit * baseinit, MPI_Comm comm) 
{
    /* plan for the variable PMs; keep in mind we do variable
     * mesh size (PM resolution). We plan them at the begining of the run
     * in theory we can use some really weird PM resolution as function
     * of time, but the parameter files / API need to support this.
     * */
    int size;
    for (size = 0; vpminit[size].pm_nc_factor > 0; size ++)
        continue;

    VPM * vpm = malloc(sizeof(VPM) * (size + 1));
    int i;
    for (i = 0; i < size; i ++) {
        vpm[i].end = 0; /* not the last item */

        vpm[i].pm_nc_factor = vpminit[i].pm_nc_factor;
        vpm[i].a_start = vpminit[i].a_start;

        PMInit pminit = *baseinit;
        pminit.Nmesh = (int)(baseinit->Nmesh * vpm[i].pm_nc_factor);
        pm_init(&vpm[i].pm, &pminit, comm);
        if(pm_unbalanced(&vpm[i].pm)) {
            fastpm_raise(-1, "PM mesh is not divided by the process mesh.\n"
                "(Nmesh[0] = %d / Nproc[0] = %d) x(Nmesh[1] = %d / Nproc[1] = %d)\n",
                "Fix this by changing the number of ranks.",
                pm_nmesh(&vpm[i].pm)[0],
                pm_nproc(&vpm[i].pm)[0],
                pm_nmesh(&vpm[i].pm)[1],
                pm_nproc(&vpm[i].pm)[1]);
        }
    }
    /* the end of the list */
    vpm[i].end = 1;
    return vpm;
}
void vpm_free (VPM * vpm) {
    int i = 0;
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
