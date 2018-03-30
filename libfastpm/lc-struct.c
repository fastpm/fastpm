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

double rad_to_degree=180./M_PI;

static int
fastpm_smesh_init_common(FastPMSMesh * mesh,
        FastPMLightCone * lc,
        size_t np_upper,
        double * a, double Na)
{
    mesh->started = 0;
    mesh->event_handlers = NULL;
    mesh->lc = lc;
    mesh->a = malloc(sizeof(double) * Na);
    mesh->z = malloc(sizeof(double) * Na);
    mesh->Na = Na;
    mesh->np_upper = np_upper;
    size_t i;

    for(i = 0; i < Na; i ++) {
        mesh->a[i] = a[i];
        mesh->z[i] = lc->speedfactor * HorizonDistance(a[i], lc->horizon);
    }

    fastpm_store_init(mesh->last.p, mesh->np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_DENSITY
            | PACK_TIDAL
            | PACK_AEMIT,
            FASTPM_MEMORY_STACK);

    mesh->last.a_f = 0;

}

void
fastpm_smesh_init_plane(FastPMSMesh * mesh,
        FastPMLightCone * lc,
        size_t np_upper,
        double (*xy)[2], size_t Nxy,
        double * a, size_t Na)
{
    mesh->type = FASTPM_SMESH_PLANE;
    mesh->Nxy = Nxy;

    fastpm_smesh_init_common(mesh, lc, np_upper, a, Na);

    mesh->xy = malloc(sizeof(double) * Nxy * 2);

    size_t i;

    for(i = 0; i < Nxy; i ++) {
        fastpm_info("init plane %ld %g %g \n",i,xy[i][0],xy[i][1]);
        mesh->xy[i][0] = xy[i][0];
        mesh->xy[i][1] = xy[i][1];
    }
}

void
fastpm_smesh_init_sphere(FastPMSMesh * mesh,
        FastPMLightCone * lc,
        size_t np_upper,
        double * ra, double * dec, size_t Nxy,
        double * a, size_t Na)
{
    mesh->type = FASTPM_SMESH_SPHERE;
    mesh->Nxy = Nxy;

    fastpm_smesh_init_common(mesh, lc, np_upper, a, Na);

    size_t i;

    mesh->vec = malloc(sizeof(double) * Nxy * 3);
    mesh->ra = malloc(sizeof(double) * Nxy);
    mesh->dec = malloc(sizeof(double) * Nxy);

    for(i = 0; i < Nxy; i ++) {
        mesh->vec[i][0] = cos(dec[i] / rad_to_degree) * cos(ra[i] / rad_to_degree);
        mesh->vec[i][1] = cos(dec[i] / rad_to_degree) * sin(ra[i] / rad_to_degree);
        mesh->vec[i][2] = sin(dec[i] / rad_to_degree);
        mesh->ra[i] = ra[i];
        mesh->dec[i] = dec[i];
    }
}

void
fastpm_smesh_destroy(FastPMSMesh * mesh)
{
    switch(mesh->type) {
        case FASTPM_SMESH_PLANE:
            free(mesh->xy);
            break;
        case FASTPM_SMESH_SPHERE:
            free(mesh->ra);
            free(mesh->dec);
            free(mesh->vec);
        break;
    }
    free(mesh->a);
    free(mesh->z);
    fastpm_destroy_event_handlers(&mesh->event_handlers);
    fastpm_store_destroy(mesh->last.p);
}

void
fastpm_smesh_select_active(FastPMSMesh * mesh,
        double a0, double a1,
        FastPMStore * q
    )
{
    size_t Na = 0;

    size_t i;
    for (i = 0; i < mesh->Na; i ++) {
        if(mesh->a[i] >= a0 && mesh->a[i] < a1) {
            Na ++;
        }
    }

    fastpm_store_init(q, mesh->np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_DENSITY
            | PACK_TIDAL
            | PACK_AEMIT,
            FASTPM_MEMORY_HEAP
    );

    if((size_t) Na * (size_t) mesh->Nxy >= q->np_upper) {
        fastpm_raise(0, "More mesh points requested than np_upper (%d * %d > %td)\n", Na, mesh->Nxy, q->np_upper);
    }

    size_t j = 0;
    size_t k = 0;
    size_t m=0;
    size_t n=0;
    double x_temp[4];
    x_temp[3]=1;
    for(j = 0; j < mesh->Nxy; j ++) {
        for(k = 0; k < mesh->Na; k ++) {
            /* cast to float because aemit is float */
            float aemit = mesh->a[k];
            if(aemit >= a0 && aemit < a1) {
                switch(mesh->type) {
                    case FASTPM_SMESH_PLANE:
                        x_temp[0]=mesh->xy[j][0];
                        x_temp[1]=mesh->xy[j][1];
                        x_temp[2]=mesh->z[k];
                        break;
                    case FASTPM_SMESH_SPHERE:
                        x_temp[0] = mesh->vec[j][0] * mesh->z[k];
                        x_temp[1] = mesh->vec[j][1] * mesh->z[k];
                        x_temp[2] = mesh->vec[j][2] * mesh->z[k];
                        break;
                }
                double xo[4];
                fastpm_gldot(mesh->lc->glmatrix_inv, x_temp, xo);
                for(m = 0; m < 3; m ++) {
                    q->x[q->np][m] = xo[m];
                }
                q->aemit[q->np] = aemit;
                q->np++;
            }
        }
    }
}

int
fastpm_smesh_compute_potential(
        FastPMSMesh * mesh,
        PM * pm,
        FastPMGravity * gravity,
        FastPMFloat * delta_k,
        double a_f,
        double a_n)
{
    double z1, z2; /* z in xyz */

    FastPMStore p_new_now[1];
    FastPMStore p_last_now[1];

    fastpm_smesh_select_active(mesh, a_f, a_n, p_new_now);

    /* create a proxy of p_last_then with the same position,
     * but new storage space for the potential variables */
    fastpm_store_init(p_last_now, mesh->np_upper,
                    mesh->last.p->attributes & ~ PACK_POS & ~ PACK_AEMIT,
                    /* skip pos, we'll use an external reference next line*/
                    FASTPM_MEMORY_HEAP
                    );
    p_last_now->np = mesh->last.p->np;
    p_last_now->x = mesh->last.p->x;
    p_last_now->aemit = mesh->last.p->aemit;
    fastpm_info("p_last_now->np = %td\n", p_last_now->np);

    FastPMFloat * canvas = pm_alloc(pm); /* Allocates memory and returns success */

    FastPMPainter reader[1];
    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);


    /*XXX Following is almost a repeat of potential calc in fastpm_gravity_calculate, though positions are different*/

    int d;
    enum FastPMPackFields ACC[] = {
                 PACK_DENSITY,
                 PACK_POTENTIAL,
                 PACK_TIDAL_XX, PACK_TIDAL_YY, PACK_TIDAL_ZZ,
                 PACK_TIDAL_XY, PACK_TIDAL_YZ, PACK_TIDAL_ZX
                };

    PMGhostData * pgd_last_now = pm_ghosts_create(pm, p_last_now, PACK_POS, NULL);
    PMGhostData * pgd_new_now = pm_ghosts_create(pm, p_new_now, PACK_POS, NULL);

    for(d = 0; d < 8; d ++) {
        CLOCK(transfer);
        gravity_apply_kernel_transfer(gravity, pm, delta_k, canvas, ACC[d]);
        LEAVE(transfer);

        CLOCK(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        CLOCK(readout);
        fastpm_readout_local(reader, canvas, p_last_now, p_last_now->np + pgd_last_now->nghosts, NULL, ACC[d]);
        fastpm_readout_local(reader, canvas, p_new_now, p_new_now->np + pgd_new_now->nghosts, NULL, ACC[d]);
        LEAVE(readout);

        CLOCK(reduce);
        pm_ghosts_reduce(pgd_last_now, ACC[d]);
        pm_ghosts_reduce(pgd_new_now, ACC[d]);
        LEAVE(reduce);

    }

    pm_ghosts_free(pgd_new_now);
    pm_ghosts_free(pgd_last_now);

    pm_free(pm, canvas);

    /* last.a_f is when the potential is last updated */
    double G_then = HorizonGrowthFactor(mesh->last.a_f, mesh->lc->horizon);
    double G_now = HorizonGrowthFactor(a_f, mesh->lc->horizon);

    /* interpolate potential and tidal field between the range, and convert the unit, into
     * last.p, which is then emitted as an event. */
    #define INTERP(field) \
        p_last_then->field = (( \
            (G_now - G_emit) * p_last_then->field \
          + (G_emit - G_then) * p_last_now->field \
            ) / (G_now - G_then))

    ptrdiff_t i;

    double potfactor = 1.5 * mesh->lc->cosmology->OmegaM / (HubbleDistance * HubbleDistance);
    FastPMStore * p_last_then = mesh->last.p;

    for(i = 0; i < p_last_now->np; i ++) {
        float a_emit = p_last_now->aemit[i];
        if(a_emit < mesh->last.a_f || a_emit >= a_f) {
            fastpm_raise(-1, " out of bounds. a_emit = %g should be between %g and %g", a_emit, mesh->last.a_f, a_f);
        }
        double G_emit = HorizonGrowthFactor(a_emit, mesh->lc->horizon);

        INTERP(potential[i]);
        p_last_then->potential[i] *= potfactor / a_emit;
        INTERP(rho[i]);
        int j;
        for(j = 0; j < 6; j ++) {
            INTERP(tidal[i][j]);
            p_last_then->tidal[i][j] *= potfactor / a_emit;
        }
    }
    /* p_last_now is no longer useful after interpolation */
    fastpm_store_destroy(p_last_now);

    /* a portion of light cone is ready between a0 and a1 */
    FastPMLCEvent lcevent[1];
    lcevent->p = mesh->last.p;
    lcevent->a0 = mesh->last.a_f;
    lcevent->a1 = a_f;
    lcevent->is_first = mesh->last.a_f == 0;

    fastpm_emit_event(mesh->event_handlers,
            FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEvent*) lcevent, mesh);

    /* last.a_f is when the potential is last updated */
    mesh->last.a_f = a_f;

    /* copy the new into the last. new is on the tack; last is on the heap. */
    fastpm_store_destroy(mesh->last.p);
    fastpm_store_init(mesh->last.p,
                    p_new_now->np_upper,
                    p_new_now->attributes,
                    FASTPM_MEMORY_STACK);

    fastpm_info("p_new_now->np = %td\n", p_new_now->np);

    fastpm_store_copy(p_new_now, mesh->last.p);


    fastpm_store_destroy(p_new_now);

    return 0;
}

