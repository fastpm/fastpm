#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/cosmology.h>
#include <fastpm/lightcone.h>
#include <fastpm/logging.h>

static Cosmology CP(FastPMSolver * fastpm) {
    Cosmology c = {
        .OmegaM = fastpm->omega_m,
        .OmegaLambda = 1 - fastpm->omega_m,
    };
    return c;
}

void
fastpm_lc_init(FastPMLightCone * lc, double speedfactor, FastPMSolver * fastpm, size_t np_upper)
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    lc->fastpm = fastpm;
    lc->p = malloc(sizeof(FastPMStore));
    /* Allocation */

    int size = 8192;
    Cosmology c = CP(fastpm);

    lc->EventHorizonTable.size = size;
    lc->EventHorizonTable.Dc = malloc(sizeof(double) * size);
    
    /* GSL init solver */
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    lc->gsl = gsl_root_fsolver_alloc(T);
    
    int i;

    for(i = 0; i < lc->EventHorizonTable.size; i ++) {
        double a = 1.0 * i / (lc->EventHorizonTable.size - 1);
        lc->EventHorizonTable.Dc[i] = speedfactor * HubbleDistance * ComovingDistance(a, c);
    }

    fastpm_store_init(lc->p);
    fastpm_store_alloc(lc->p, np_upper, PACK_ID | PACK_POS | PACK_VEL | PACK_AEMIT);
}

void
fastpm_lc_destroy(FastPMLightCone * lc)
{
    /* Free */
    fastpm_store_destroy(lc->p);
    free(lc->EventHorizonTable.Dc);
    free(lc->p);
    /* GSL destroy solver */
    gsl_root_fsolver_free(lc->gsl);
}

double
fastpm_lc_horizon(FastPMLightCone * lc, double a)
{
    /* It may be worth to switch to log_a interpolation, but it only matters
     * at very high z (~ z = 9). */

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

struct funct_params {
    FastPMLightCone *lc;
    FastPMStore * p;
    FastPMDrift * drift;
    ptrdiff_t i;
};

static double
funct(double a, void *params)
{
    struct funct_params *Fp = (struct funct_params *) params;

    FastPMLightCone *lc = Fp->lc;
    FastPMStore * p = Fp->p;
    FastPMDrift * drift = Fp->drift;
    ptrdiff_t i = Fp->i;
    double xo[3];

    fastpm_drift_one(drift, p, i, xo, a);

    /*XXX: may need to worry about periodic boundary */
    return xo[2] - fastpm_lc_horizon(lc, a);
}

static int
_fastpm_lc_intersect_one(FastPMLightCone * lc,
        struct funct_params * params,
        ptrdiff_t i,
        double * solution)
{
    params->i = i;

    int status;
    int iter = 0, max_iter;
    double r, x_lo, x_hi, eps;

    /* Reorganize to struct later */
    x_lo = params->drift->ai;
    x_hi = params->drift->af;
    max_iter = 100;
    eps = 1e-7;

    gsl_function F;

    F.function = &funct;
    F.params = params;

    status = gsl_root_fsolver_set(lc->gsl, &F, x_lo, x_hi);

    if(status == GSL_EINVAL || status == GSL_EDOM) { 
        /** Error in value or out of range **/
        return 0; 
    } 

    do
    {
        iter++;
        //
        // Debug printout #1
        //if(iter == 1) {
        //fastpm_info("ID | [x_lo, x_hi] | r | funct(r) | x_hi - x_lo\n");
        //}
        //

        status = gsl_root_fsolver_iterate(lc->gsl);
        r = gsl_root_fsolver_root(lc->gsl);

        x_lo = gsl_root_fsolver_x_lower(lc->gsl);
        x_hi = gsl_root_fsolver_x_upper(lc->gsl);

        status = gsl_root_test_interval(x_lo, x_hi, eps, 0.0);
        //
        //Debug printout #2
        //fastpm_info("%5d [%.7f, %.7f] %.7f %.7f %.7f\n", iter, x_lo, x_hi, r, funct(r, &params), x_hi - x_lo);
        //

        if(status == GSL_SUCCESS) {
            *solution = r;
            //
            // Debug printout #3.1
            //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 1);
            //
            return 1;
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    //
    // Debug printout #3.2
    //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 0);
    //

    return 0;
}

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMKick * kick, FastPMStore * p)
{
    struct funct_params params = {
        .lc = lc, 
        .drift = drift,
        .p = p,
        .i = 0,
    };
    ptrdiff_t i;

    for(i = 0; i < p->np; i ++) {
        double a_emit = 0;
        if(0 == _fastpm_lc_intersect_one(lc, &params, i, &a_emit)) continue;
        /* A solution is found */
        /* move the particle and store it. */
        ptrdiff_t next = lc->p->np;
        if(next == lc->p->np_upper) {
            fastpm_raise(-1, "Too many particles in the light cone");
        }
        double xo[3];
        float vo[3];
        fastpm_drift_one(drift, p, i, xo, a_emit);
        fastpm_kick_one(kick, p, i, vo, a_emit);
        int d;
        for(d = 0; d < 3; d ++) {
            lc->p->x[next][d] = xo[d];
            /* XXX: convert units? */
            lc->p->v[next][d] = vo[d];
        }
        lc->p->id[next] = p->id[i];
        lc->p->aemit[next] = a_emit;
        lc->p->np ++;
    }
    return 0;
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
