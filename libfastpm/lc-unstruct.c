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
#include <gsl/gsl_linalg.h>

#include "pmpfft.h"
#include "pmghosts.h"

static void
_matrix_invert(double (*matrix)[4],
                    double (*matrix_inv)[4],
                    int matrix_size) 
                                                                        //square matrix.. 
                                                                        //matrix_size fixed to 4 
{
    double temp[16];

    memcpy(temp, matrix, sizeof(double) * 16);

    gsl_matrix_view m = gsl_matrix_view_array(temp, matrix_size, matrix_size);
    gsl_matrix_view inv = gsl_matrix_view_array(&matrix_inv[0][0],matrix_size, matrix_size);
    gsl_permutation * p = gsl_permutation_alloc (matrix_size);

    int signum;
    gsl_linalg_LU_decomp (&m.matrix, p, &signum);    
    gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);

    gsl_permutation_free (p);
}

void
fastpm_lc_init(FastPMLightCone * lc)
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    /* Allocation */
    lc->horizon = malloc(sizeof(FastPMHorizon));
    fastpm_horizon_init(lc->horizon, lc->cosmology);
    _matrix_invert(lc->glmatrix,lc->glmatrix_inv,4);

}

void
fastpm_lc_destroy(FastPMLightCone * lc)
{
    fastpm_horizon_destroy(lc->horizon);
    free(lc->horizon);
}

void
fastpm_usmesh_init(FastPMUSMesh * mesh, FastPMLightCone * lc,
            size_t np_upper,
            double (*tileshifts)[3],
            int ntiles,
            double amin,
            double amax)
{

    mesh->amin = amin;
    mesh->amax = amax;
    mesh->lc = lc;
    mesh->tileshifts = malloc(sizeof(tileshifts[0]) * ntiles);
    mesh->ntiles = ntiles;
    mesh->is_first = 1;

    memcpy(mesh->tileshifts, tileshifts, sizeof(tileshifts[0]) * ntiles);

    mesh->event_handlers = NULL;
    mesh->p = malloc(sizeof(FastPMStore));
    /* for saving the density with particles */
    fastpm_store_init(mesh->p, np_upper,
                  PACK_ID | PACK_POS | PACK_VEL
                | PACK_AEMIT,
                FASTPM_MEMORY_HEAP
    );
}

void fastpm_usmesh_destroy(FastPMUSMesh * mesh)
{
    fastpm_destroy_event_handlers(&mesh->event_handlers);
    fastpm_store_destroy(mesh->p);
    free(mesh->tileshifts);
    free(mesh->p);
}

struct funct_params {
    FastPMLightCone *lc;
    FastPMStore * p;
    FastPMDriftFactor * drift;
    ptrdiff_t i;
    double tileshift[4];
    double a1;
    double a2;
    void * context;
};

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
    fastpm_gldot(lc->glmatrix, xi, xo);

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

    return distance - lc->speedfactor * HorizonDistance(a, lc->horizon);
}

static int
_fastpm_usmesh_intersect_one(FastPMUSMesh * mesh,
        struct funct_params * params,
        ptrdiff_t i,
        double * solution)
{
    params->i = i;

    return fastpm_horizon_solve(mesh->lc->horizon,
        params->context,
        solution,
        params->a1, params->a2,
        funct, params);
}

static double
zangle(double * x) {
    double dxy = 0;
    double dz = x[2];
    dxy = x[0] * x[0] + x[1] * x[1];

    double rt = atan2(sqrt(dxy), dz) / M_PI * 180.;
    if (rt < 0) rt += 360.;
    return rt;
}

static int
_in_octant(int i, double vec[3], double tol)
{
    static const int signs[][3] = {
        {1, 1, 1},
        {1, 1, -1},
        {1, -1, 1},
        {1, -1, -1},
        {-1, 1, 1},
        {-1, 1, -1},
        {-1, -1, 1},
        {-1, -1, -1},
    };
    int d;
    for(d = 0; d < 3; d ++) {
        double s = vec[d] * signs[i][d];
        if(s < -tol) {
            return 0;
        }
    }
    return 1;
}
int
fastpm_lc_inside(FastPMLightCone * lc, double vec[3])
{
    if(lc->fov > 0) {
        int r;
        if(lc->fov < 360) {
            r = zangle(vec) <= lc->fov * 0.5;
        } else {
            r = 1;
        }
        if(r) {
            double norm = 0;
            int d;
            for(d = 0; d < 3; d ++) {
                norm += vec[d] * vec[d];
            }
            norm = sqrt(norm);
            int i;
            for(i = 0; i < 8; i ++) {
                if(lc->octants[i] && _in_octant(i, vec, lc->tol * norm)) return 1;
            }
            return 0;
        } else {
            return 0;
        }
    } else {
        /* smesh takes care of this */
        return 1;
    }
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
fastpm_usmesh_intersect_tile(FastPMUSMesh * mesh, double * tileshift,
        double a1,
        double a2,
        FastPMDriftFactor * drift,
        FastPMKickFactor * kick,
        FastPMStore * p,
        FastPMStore * pout
)
{
    FastPMLightCone * lc = mesh->lc;
    struct funct_params params = {
        .lc = lc,
        .drift = drift,
        .p = p,
        .i = 0,
        .a1 = a1, 
        .a2 = a2,
    };

    int a1_is_outside = (params.a1 > mesh->amax) || (params.a1 < mesh->amin);
    int a2_is_outside = (params.a2 > mesh->amax) || (params.a2 < mesh->amin);

    if(a1_is_outside && a2_is_outside) {
        return 0;
    }

    {
        int d;
        for(d = 0; d < 3; d ++) {
            params.tileshift[d] = tileshift[d];
        }
    }

    fastpm_info("Considering Tile %g %g %g\n",
        params.tileshift[0],
        params.tileshift[1],
        params.tileshift[2]);

    params.tileshift[3] = 0;

    ptrdiff_t i;

    #pragma omp parallel firstprivate(params)
    {
        params.context = fastpm_horizon_solve_start();
        #pragma omp for
        for(i = 0; i < p->np; i ++) {
            double a_emit = 0;
            if(0 == _fastpm_usmesh_intersect_one(mesh, &params, i, &a_emit)) continue;
            
            /* the event is outside the region we care, skip */
            if(a_emit > mesh->amax || a_emit < mesh->amin) continue;

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
            fastpm_gldot(lc->glmatrix, xi, xo);

            /* does it fall into the field of view? */
            if(!fastpm_lc_inside(lc, xo)) continue;

            /* A solution is found */
            /* move the particle and store it. */
            ptrdiff_t next;
            #pragma omp atomic capture
                next = pout->np++;

            if(next == pout->np_upper) {
                fastpm_raise(-1, "Too many particles in the light cone; limit = %td\n", pout->np_upper);
            }

            /* copy the position if desired */
            if(pout->x) {
                for(d = 0; d < 3; d ++) {
                    pout->x[next][d] = xo[d];
                }
            }

            float vo[4];
            float vi[4];
            if(p->v) {
                /* can we kick? if we are using a fixed grid there is no v */
                fastpm_kick_one(kick, p, i, vi, a_emit);
                vi[3] = 0;
                /* transform the coordinate */
                fastpm_gldotf(lc->glmatrix, vi, vo);

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
        }
        fastpm_horizon_solve_end(params.context);
    }

    return 0;
}

int
fastpm_usmesh_intersect(FastPMUSMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm)
{
    double a1 = drift->ai > drift->af ? drift->af: drift->ai;
    double a2 = drift->ai > drift->af ? drift->ai: drift->af;

    /* for each tile */
    int t;
    for(t = 0; t < mesh->ntiles; t ++) {
        fastpm_usmesh_intersect_tile(mesh, &mesh->tileshifts[t][0],
                a1, a2,
                drift, kick,
                fastpm->p,
                mesh->p); /*Store particle to get density*/

    }

    /* a portion of light cone is ready between a0 and a1 */
    FastPMLCEvent lcevent[1];

    lcevent->p = mesh->p;
    lcevent->a0 = a1;
    lcevent->a1 = a2;
    lcevent->is_first = mesh->is_first;

    fastpm_emit_event(mesh->event_handlers,
            FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEvent*) lcevent, mesh);

    /* now purge the store. */
    mesh->p->np = 0;
    mesh->is_first = 0;

    return 0;
}
