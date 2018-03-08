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
fastpm_struct_mesh_init_plane(FastPMStructuredMesh * mesh,
        double (*xy)[2], size_t Nxy,
        double * z, size_t Nz)
{
    mesh->type = FASTPM_STRUCT_MESH_PLANE;
    mesh->z = malloc(sizeof(double) * Nz);
    mesh->Nz = Nz;
    size_t i;

    for(i = 0; i < Nz; i ++) {
        mesh->z[i] = z[i];
    }

    mesh->xy = malloc(sizeof(double) * Nxy * 2);
    for(i = 0; i < Nxy; i ++) {
        mesh->xy[i][0] = xy[i][0];
        mesh->xy[i][1] = xy[i][1];
    }
}

void
fastpm_struct_mesh_init_sphere(FastPMStructuredMesh * mesh,
        double * ra, double * dec, size_t Npix,
        double * z, size_t Nz)
{
    mesh->type = FASTPM_STRUCT_MESH_SPHERE;
    mesh->z = malloc(sizeof(double) * Nz);
    mesh->Nz = Nz;
    size_t i;

    for(i = 0; i < Nz; i ++) {
        mesh->z[i] = z[i];
    }

    mesh->vec = malloc(sizeof(double) * Npix * 3);
    mesh->ra = malloc(sizeof(double) * Npix);
    mesh->dec = malloc(sizeof(double) * Npix);

    for(i = 0; i < Npix; i ++) {
        /* FIXME: does this look correct? */
        mesh->vec[i][0] = cos(dec[i]) * cos(ra[i]);
        mesh->vec[i][1] = cos(dec[i]) * sin(ra[i]);
        mesh->vec[i][2] = sin(dec[i]);
        mesh->ra[i] = ra[i];
        mesh->dec[i] = dec[i];
    }
}

void
fastpm_struct_mesh_select_active(FastPMStructuredMesh * mesh,
        double z0, double z1,
        FastPMStore * q
    )
{
    size_t Nz = 0;

    size_t i;
    for (i = 0; i < mesh->Nz; i ++) {
        if(mesh->z[i] >= z0 && mesh->z[i] < z1) {
            Nz ++;
        }
    }

    /* 2 is fudge factor for ghosts */
    size_t np_upper = mesh->Nxy * Nz * 2;

    fastpm_store_init(q, np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_TIDAL
            | PACK_AEMIT,
            FASTPM_MEMORY_STACK /* from the stack so it won't interfere with the last. */
    );
    size_t j = 0;
    size_t k = 0;
    for(j = 0; j < mesh->Nxy; j ++) {
        for(k = 0; k < mesh->Nz; k ++) {
            if(mesh->z[i] >= z0 && mesh->z[i] < z1) {
                switch(mesh->type) {
                    case FASTPM_STRUCT_MESH_PLANE:
                        q->x[q->np][0] = mesh->xy[j][0];
                        q->x[q->np][1] = mesh->xy[j][1];
                        q->x[q->np][2] = mesh->z[j];
                        break;
                    case FASTPM_STRUCT_MESH_SPHERE:
                        q->x[q->np][0] = mesh->vec[j][0] * mesh->z[j];
                        q->x[q->np][1] = mesh->vec[j][1] * mesh->z[j];
                        q->x[q->np][2] = mesh->vec[j][2] * mesh->z[j];
                        break;
                }
                q->np++;
            }
        }
    }
}

void
fastpm_struct_mesh_intersect(FastPMSolver * fastpm, FastPMTransitionEvent * event, FastPMStructuredMesh * mesh)
{


}

int
fastpm_struct_mesh_compute_potential(FastPMSolver * fastpm,
        FastPMForceEvent * event,
        FastPMStructuredMesh * mesh)
{
    PM * pm = fastpm->pm;

    FastPMFloat * delta_k = event->delta_k;
    FastPMFloat * canvas = pm_alloc(pm); /* Allocates memory and returns success */
    FastPMGravity * gravity = fastpm->gravity;
    FastPMPainter reader[1];
    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);

    double z0, z1; /* z in xyz */

    z0 = mesh->last.z1;
    z1 = HorizonDistance(event->a_f, mesh->lc->horizon);

    FastPMStore p_new_now[1];
    FastPMStore p_last_now[1];

    fastpm_struct_mesh_select_active(mesh, z0, z1, p_new_now);

    /* create a proxy of p_last_then with the same position,
     * but new storage space for the potential variables */
    fastpm_store_init(p_last_now, mesh->last.p->np_upper,
                    mesh->last.p->attributes & ~ PACK_POS, /* skip pos, we'll use an external reference next line*/
                    FASTPM_MEMORY_STACK
                    );
    p_last_now->x = mesh->last.p->x;

    /*XXX Following is almost a repeat of potential calc in fastpm_gravity_calculate, though positions are different*/

    int d;
    int ACC[] = {
                 PACK_POTENTIAL,
                 PACK_TIDAL_XX, PACK_TIDAL_YY, PACK_TIDAL_ZZ,
                 PACK_TIDAL_XY, PACK_TIDAL_YZ, PACK_TIDAL_ZX
                };

    PMGhostData * pgd_last_now = pm_ghosts_create(pm, p_last_now, PACK_POS, NULL);
    PMGhostData * pgd_new_now = pm_ghosts_create(pm, p_new_now, PACK_POS, NULL);

    for(d = 0; d < 7; d ++) {
        CLOCK(transfer);
        gravity_apply_kernel_transfer(gravity, pm, delta_k, canvas, d+3);
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

    pm_free(pm, canvas);

    pm_ghosts_free(pgd_new_now);
    pm_ghosts_free(pgd_last_now);

    /* FIXME: read off the potential, probably emit another event to write out
     * p_last_now? */

    FastPMLCEvent lcevent[1];
    lcevent->p = mesh->last.p;
    lcevent->z0 = mesh->last.z0;
    lcevent->z1 = mesh->last.z1;

    fastpm_emit_event(mesh->event_handlers,
            FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
            (FastPMEvent*) lcevent, mesh);

    fastpm_store_destroy(p_last_now);
    fastpm_store_destroy(mesh->last.p);

    /* copy the new into the last. new is on the tack; last is on the heap. */
    fastpm_store_init(mesh->last.p,
                    p_new_now->np_upper,
                    p_new_now->attributes,
                    FASTPM_MEMORY_HEAP);

    fastpm_store_copy(p_new_now, mesh->last.p);

    mesh->last.z1 = z1;
    mesh->last.z0 = z0;

    fastpm_store_destroy(p_new_now);

    return 0;
}

