#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

#include "libfastpm.h"
#include "pmpfft.h"
#include "pmstore.h"
#include "pmghosts.h"

void
fastpm_paint(PM * pm, PMStore * p, FastPMFloat * delta_x, FastPMFloat * delta_k) 
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
fastpm_dump(PM * pm , char * filename, FastPMFloat *data) 
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
