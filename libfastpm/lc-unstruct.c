#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/lc-unstruct.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmghosts.h"

void
fastpm_lc_init(FastPMLightCone * lc, FastPMSolver * fastpm,
                double (*tileshifts)[3], int ntiles
                )
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    lc->unstruct = malloc(sizeof(FastPMStore));

    /* Allocation */

    lc->tileshifts = malloc(sizeof(tileshifts[0]) * ntiles);
    lc->ntiles = ntiles;

    memcpy(lc->tileshifts, tileshifts, sizeof(tileshifts[0]) * ntiles);

    lc->horizon = malloc(sizeof(FastPMHorizon));
    fastpm_horizon_init(lc->horizon, lc->cosmology);

    /* for saving the density with particles */
    fastpm_store_init(lc->unstruct, fastpm->p->np_upper,
                  PACK_ID | PACK_POS | PACK_VEL
                | PACK_AEMIT
                | (fastpm->p->potential?PACK_POTENTIAL:0)
    );

}

void
fastpm_lc_destroy(FastPMLightCone * lc)
{
    /* Free */
    fastpm_store_destroy(lc->unstruct);

    free(lc->tileshifts);
    free(lc->unstruct);

    fastpm_horizon_destroy(lc->horizon);
    free(lc->horizon);
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

    xi[3] = 1;
    if(p->v) {
        fastpm_drift_one(drift, p, i, xi, a);
    } else {
        for(d = 0; d < 3; d ++) {
            xi[d] = p->x[i][d];
        }
    }
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

    return distance - HorizonDistance(a, lc->horizon);
}

static int
_fastpm_lc_intersect_one(FastPMLightCone * lc,
        struct funct_params * params,
        ptrdiff_t i,
        double * solution)
{
    params->i = i;

    return fastpm_horizon_solve(lc->horizon,
        solution,
        params->a1, params->a2,
        funct, params);
}

static double
zangle(double * x) {
    double dxy = 0;
    double dz = x[2];
    dxy = x[0] * x[0] + x[1] * x[1];

    return atan2(sqrt(dxy), dz) / M_PI * 180.;
}

/* FIXME:
 * the function shall take ai, af as input,
 *
 * the function shall be able to interpolate potential as
 * well as position and velocity.
 * We need a more general representation of '*drift' and '*kick'.
 *
 * */
static int
fastpm_lc_intersect_tile(FastPMLightCone * lc, int tile,
        FastPMDriftFactor * drift,
        FastPMKickFactor * kick,
        FastPMStore * p,
        FastPMStore * pout
)
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
        ptrdiff_t next = pout->np;
        if(next == pout->np_upper) {
            fastpm_raise(-1, "Too many particles in the light cone");
        }
        double xi[4];
        double xo[4];
        int d;

        xi[3] = 1;
        if(p->v) {
            /* can we drift? if we are using a fixed grid there is no v. */
            fastpm_drift_one(drift, p, i, xi, a_emit);
        } else {
            for(d = 0; d < 3; d ++) {
                xi[d] = p->x[i][d];
            }
        }
        for(d = 0; d < 4; d ++) {
            xi[d] += params.tileshift[d];
        }
        /* transform the coordinate */
        gldot(lc->glmatrix, xi, xo);

        /* does it fall into the field of view? */
        if(lc->fov > 0 && zangle(xo) > lc->fov * 0.5) continue;

        /* copy the position if desired */
        if(pout->x) {
            for(d = 0; d < 3; d ++) {
                pout->x[next][d] = xo[d];
            }
        }

        float vo[3];
        float vi[3];
        if(p->v) {
            /* can we kick? if we are using a fixed grid there is no v */
            fastpm_kick_one(kick, p, i, vi, a_emit);
            /* transform the coordinate */
            gldotv(lc->glmatrix, vi, vo);

            if(pout->v) {
                for(d = 0; d < 3; d ++) {
                    /* convert to peculiar velocity a dx / dt in kms */
                    pout->v[next][d] = vo[d] * HubbleConstant / a_emit;
                }
            }
        }
        if(pout->id)
            pout->id[next] = p->id[i];
        if(pout->aemit)
            pout->aemit[next] = a_emit;

        double potfactor = 1.5 * lc->cosmology->OmegaM / (HubbleDistance * HubbleDistance);
        /* convert to dimensionless potential */
        if(pout->potential)
            pout->potential[next] = p->potential[i] / a_emit * potfactor;

        if(pout->tidal) {
            for(d = 0; d < 6; d++) {
                pout->tidal[next][d] = p->tidal[i][d] / a_emit * potfactor;
            }
        }
        pout->np ++;
    }
    return 0;
}

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm)
{
    /* for each tile */
    int t;
    for(t = 0; t < lc->ntiles; t ++) {
        fastpm_lc_intersect_tile(lc, t, drift, kick, fastpm->p, lc->unstruct);/*Store particle to get density*/

    }
    return 0;
}
