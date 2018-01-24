#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/lightcone.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmghosts.h"

void
fastpm_lc_init(FastPMLightCone * lc, FastPMSolver * fastpm,
                double (*tileshifts)[3], int ntiles
                )
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    lc->p = malloc(sizeof(FastPMStore));
    lc->q = malloc(sizeof(FastPMStore));
    lc->p0 = malloc(sizeof(FastPMStore));

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

    if(lc->compute_potential) {
        /* p0 is the lagrangian _position */
        fastpm_store_init(lc->p0, fastpm->p->np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_TIDAL
        );

        fastpm_store_set_lagrangian_position(lc->p0, fastpm->basepm, NULL, NULL);

        /* for saving the lagrangian sampling of potential */
        fastpm_store_init(lc->q, fastpm->p->np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_TIDAL
            | PACK_AEMIT);
    }
    /* for saving the density with particles */
    fastpm_store_init(lc->p, fastpm->p->np_upper,
                  PACK_ID | PACK_POS | PACK_VEL
                | PACK_AEMIT
                | (fastpm->p->potential?PACK_POTENTIAL:0)
    );

}

void
fastpm_lc_destroy(FastPMLightCone * lc)
{
    /* Free */
    fastpm_store_destroy(lc->p);

    if(lc->compute_potential) {
        fastpm_store_destroy(lc->q);
        fastpm_store_destroy(lc->p0);
    }
    free(lc->EventHorizonTable.Dc);
    free(lc->tileshifts);
    free(lc->p0);
    free(lc->q);
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

    return distance - fastpm_lc_horizon(lc, a);
}

int
fastpm_lc_compute_potential(FastPMSolver * fastpm,
        FastPMForceEvent * event,
        FastPMLightCone * lc)
{

    PM * pm = fastpm->pm;
    FastPMFloat * delta_k = event->delta_k;
    FastPMFloat * canvas = pm_alloc(pm);/*Allocates memory and returns success*/
    FastPMGravity * gravity = fastpm->gravity;
    FastPMPainter reader[1];
    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);

    FastPMStore * p = lc->p0;

    int d;
    int ACC[] = {
                 PACK_POTENTIAL,
                 PACK_TIDAL_XX, PACK_TIDAL_YY, PACK_TIDAL_ZZ,
                 PACK_TIDAL_XY, PACK_TIDAL_YZ, PACK_TIDAL_ZX
                };

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    for(d = 0; d < 7; d ++) {
        CLOCK(transfer);
        gravity_apply_kernel_transfer(gravity, pm, delta_k, canvas, d);
        LEAVE(transfer);

        CLOCK(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        CLOCK(readout);
        fastpm_readout_local(reader, canvas, p, p->np + pgd->nghosts, NULL, ACC[d]);
        LEAVE(readout);

        CLOCK(reduce);
        pm_ghosts_reduce(pgd, ACC[d]);
        LEAVE(reduce);

    }
    pm_free(pm, canvas);

    pm_ghosts_free(pgd);

    return 0;
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
        fastpm_lc_intersect_tile(lc, t, drift, kick, fastpm->p, lc->p);/*Store particle to get density*/

        if(lc->compute_potential)
            fastpm_lc_intersect_tile(lc, t, drift, kick, lc->p0, lc->q);/*store potential on fixed grid*/
    }
    return 0;
}
