#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/lightcone.h>
#include <fastpm/logging.h>

void
fastpm_lc_init(FastPMLightCone * lc, FastPMStore * p,
                double (*tileshifts)[3], int ntiles
                )
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    lc->p = malloc(sizeof(FastPMStore));
    /* Allocation */

    int size = 8192;

    lc->EventHorizonTable.size = size;
    lc->EventHorizonTable.Dc = malloc(sizeof(double) * size);
    lc->tileshifts = malloc(sizeof(tileshifts[0]) * ntiles);
    lc->ntiles = ntiles;

    memcpy(lc->tileshifts, tileshifts, sizeof(tileshifts[0]) * ntiles);

    /* GSL init solver */
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    lc->gsl = gsl_root_fsolver_alloc(T);

    int i;

    for(i = 0; i < lc->EventHorizonTable.size; i ++) {
        double a = 1.0 * i / (lc->EventHorizonTable.size - 1);
        lc->EventHorizonTable.Dc[i] = lc->speedfactor * HubbleDistance * ComovingDistance(a, lc->cosmology);
    }

    fastpm_store_init(lc->p, p->np_upper, PACK_ID | PACK_POS | PACK_VEL | PACK_AEMIT | (p->potential?PACK_POTENTIAL:0));
}

void
fastpm_lc_destroy(FastPMLightCone * lc)
{
    /* Free */
    fastpm_store_destroy(lc->p);
    free(lc->EventHorizonTable.Dc);
    free(lc->tileshifts);
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
    FastPMDriftFactor * drift;
    ptrdiff_t i;
    double tileshift[4];
    double a1;
    double a2;
};

static void
gldot(double glmatrix[4][4], double xi[4], double xo[4])
{
    int i, j;
    for(i = 0; i < 4; i ++) {
        xo[i] = 0;
        for(j = 0; j < 4; j ++) {
            xo[i] += glmatrix[i][j] * xi[j];
        }
    }
}

static void
gldotv(double glmatrix[4][4], float vi[3], float vo[3])
{
    int i, j;
    for(i = 0; i < 3; i ++) {
        vo[i] = 0;
        for(j = 0; j < 3; j ++) {
            vo[i] += glmatrix[i][j] * vi[j];
        }
    }
}

static double
funct(double a, void *params)
{
    struct funct_params *Fp = (struct funct_params *) params;

    FastPMLightCone *lc = Fp->lc;
    FastPMStore * p = Fp->p;
    FastPMDriftFactor * drift = Fp->drift;
    ptrdiff_t i = Fp->i;
    int d;
    double xi[4];
    double xo[4];

    fastpm_drift_one(drift, p, i, xi, a);
    xi[3] = 1;
    for(d = 0; d < 4; d ++) {
        xi[d] += Fp->tileshift[d];
    }
    /* transform the coordinate */
    gldot(lc->glmatrix, xi, xo);

    /* XXX: may need to worry about periodic boundary */
    double distance;
    if (lc->fov <= 0) {
        distance = xo[2];
    } else {
        distance = 0;
        for (d = 0; d < 3; d ++) {
            distance += xo[d] * xo[d];
        }
        distance = sqrt(distance);
    }

    return distance - fastpm_lc_horizon(lc, a);
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
    x_lo = params->a1;
    x_hi = params->a2;
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

static int
fastpm_lc_intersect_tile(FastPMLightCone * lc, int tile, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMStore * p)
{
    struct funct_params params = {
        .lc = lc, 
        .drift = drift,
        .p = p,
        .i = 0,
        .a1 = drift->ai > drift->af ? drift->af: drift->ai,
        .a2 = drift->ai > drift->af ? drift->ai: drift->af,
    };
    int d;

    for(d = 0; d < 3; d ++) {
        params.tileshift[d] = lc->tileshifts[tile][d];
    }
    fastpm_info("tileshift = %g %g %g\n",
        params.tileshift[0],
        params.tileshift[1],
        params.tileshift[2]);

    params.tileshift[3] = 0;

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
        double xi[4];
        double xo[4];
        float vo[3];
        float vi[3];
        fastpm_drift_one(drift, p, i, xi, a_emit);
        fastpm_kick_one(kick, p, i, vi, a_emit);

        int d;
        xi[3] = 1;
        for(d = 0; d < 4; d ++) {
            xi[d] += params.tileshift[d];
        }
        /* transform the coordinate */
        gldot(lc->glmatrix, xi, xo);
        /* transform the coordinate */
        gldotv(lc->glmatrix, vi, vo);

        for(d = 0; d < 3; d ++) {
            lc->p->x[next][d] = xo[d];
            /* convert to peculiar velocity a dx / dt in kms */
            lc->p->v[next][d] = vo[d] * HubbleConstant / a_emit;
        }
        lc->p->id[next] = p->id[i];
        lc->p->aemit[next] = a_emit;
        double potfactor = 1.5 * lc->cosmology->OmegaM / (HubbleDistance * HubbleDistance);
        /* convert to dimensionless potential */
        if(lc->p->potential)
            lc->p->potential[next] = p->potential[i] / a_emit * potfactor;
        lc->p->np ++;
    }
    return 0;
}

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMStore * p)
{
    /* for each tile */
    int t;
    for(t = 0; t < lc->ntiles; t ++) {
        fastpm_lc_intersect_tile(lc, t, drift, kick, p);
    }
    return 0;
}

