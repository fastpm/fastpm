#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#include <fastpm/libfastpm.h>

#include <fastpm/cosmology.h>
#include <fastpm/lightcone.h>

void
fastpm_lc_init(FastPMLightCone * lc, Cosmology CP, size_t np_upper)
{
    /* allocation */
    int size = 8192;
    lc->EventHorizonTable.size = size;
    lc->EventHorizonTable.Dc = malloc(sizeof(double) * size);
    int i;
    for(i = 0; i < lc->EventHorizonTable.size; i ++) {
        double a = 1.0 * i / (lc->EventHorizonTable.size - 1);
        lc->EventHorizonTable.Dc[i] = ComovingDistance(a, CP);
    }
}

void fastpm_lc_destroy(FastPMLightCone * lc)
{

    /* free */
    free(lc->EventHorizonTable.Dc);
}

double
fastpm_lc_horizon(FastPMLightCone * lc, double a)
{
    /* It may be worth it to switch to loga interpolation but it only matters
     * at very high z. (~ z = 9) */

    double x = a * (lc->EventHorizonTable.size - 1);
    int l = floor(x);
    int r = l + 1;
    if(r >= lc->EventHorizonTable.size) {
        return lc->EventHorizonTable.Dc[lc->EventHorizonTable.size - 1];
    }
    if(l <= 0) {
        return lc->EventHorizonTable.Dc[0];
    }
    return lc->EventHorizonTable.Dc[l] * (r - x)
         + lc->EventHorizonTable.Dc[r] * (x - l);
}

#if 0

void fastpm_lc_write(FastPMLightCone * lc, const char * filename)
{

}

double _unknown(double a, void * params)
{

    double xo[3];
    double horizon;

    horizon = fastpm_lc_horizon(lc, a);

    fastpm_drift_one(drift, p, i, xo, a);

    return horizon - xo[2];
}

int fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p, int i, double * solution)
{

    /* XXX: culling */

    gsl_fsolve(solution);

    if( *solution >= drift->a_i
     && *solution <= drift->a_f)
    {
        return 1;
    }
    return 0;
}

int fastpm_lc_scan(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p)
{
    for(i in p) {
        status = fastpm_lc_intersect(lc, drift, p, i, solution)
        if(dint is good)
            add(i, lc);
    }
}
#endif
