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

static void
fill_a(FastPMSMesh * mesh, double * a, int Na,
       double zmin, double zmax);

void
fastpm_smesh_init(FastPMSMesh * mesh, FastPMLightCone * lc, size_t np_upper, double smoothing)
{
    mesh->event_handlers = NULL;
    mesh->lc = lc;
    mesh->np_upper = np_upper;
    mesh->layers = NULL;
    mesh->smoothing = smoothing;
    fastpm_store_init(mesh->last.p, 0,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_DENSITY
            | PACK_TIDAL
            | PACK_AEMIT,
            FASTPM_MEMORY_STACK);

    mesh->last.a_f = 0;
    mesh->started = 0;
}

struct FastPMSMeshLayer *
fastpm_smesh_add_layer_common(FastPMSMesh * mesh,
        double * a, double Na)
{
    struct FastPMSMeshLayer * h;
    struct FastPMSMeshLayer * layer = malloc(sizeof(struct FastPMSMeshLayer));

    h = mesh->layers;
    mesh->layers = layer;
    layer->next = h;

    layer->a = malloc(sizeof(double) * Na);
    layer->z = malloc(sizeof(double) * Na);
    layer->Na = Na;
    size_t i;

    for(i = 0; i < Na; i ++) {
        layer->a[i] = a[i];
        layer->z[i] = mesh->lc->speedfactor * HorizonDistance(a[i], mesh->lc->horizon);
    }

    return layer;
}

void
fastpm_smesh_add_layer_plane(FastPMSMesh * mesh,
        double (*xy)[2], size_t Nxy,
        double * a, size_t Na)
{
    struct FastPMSMeshLayer * layer = fastpm_smesh_add_layer_common(mesh, a, Na);

    layer->type = FASTPM_SMESH_PLANE;
    layer->Nxy = Nxy;

    layer->xy = malloc(sizeof(double) * Nxy * 2);

    size_t i;

    for(i = 0; i < Nxy; i ++) {
        layer->xy[i][0] = xy[i][0];
        layer->xy[i][1] = xy[i][1];
    }
}

void
fastpm_smesh_add_layer_pm(FastPMSMesh * mesh,
        PM * pm, double * shift, ptrdiff_t * Nc, double amin, double amax)
{
    double * BoxSize = pm_boxsize(pm);
    /* creat a mesh with a uniform xy grid, respecting domain given by pm.
     * use a subsample ratio. 
     * (every subsample grid points) */
    if(Nc == NULL) {
        Nc = pm_nmesh(pm);
    }
    double noshift[] = {0, 0, 0};
    if(shift == NULL) {
        shift = noshift;
    }
    size_t Nxy;
    int start[2];
    int end[2];
    int d;
    int i, j;

    Nxy = 1;

    for(d = 0; d < 2; d++) {
        start[d] = pm->IRegion.start[d] * Nc[d] / pm->Nmesh[d];
        end[d] = (pm->IRegion.start[d] + pm->IRegion.size[d]) * Nc[d] / pm->Nmesh[d];
        Nxy *= end[d] - start[d];
    }

    double zmin = mesh->lc->speedfactor * HorizonDistance(amax, mesh->lc->horizon);
    double zmax = mesh->lc->speedfactor * HorizonDistance(amin, mesh->lc->horizon);

    int Na = ceil(Nc[2] / BoxSize[2] * (zmax - zmin)) + 1;

    if (Na < 0) Na = 1;

    double (*xy)[2] = malloc(sizeof(double) * 2 * Nxy);
    double (*a) = malloc(sizeof(double) * Na);

    ptrdiff_t ptr = 0;
    for(i = start[0] ; i < end[0]; i ++) {
        for(j = start[1] ; j < end[1]; j ++) {
            /* remember the observer shall be looking at the center of the mesh. */
            xy[ptr][0] = (i + shift[0]) * BoxSize[0] / Nc[0] - 0.5 * BoxSize[0];
            xy[ptr][1] = (j + shift[0]) * BoxSize[1] / Nc[1] - 0.5 * BoxSize[1];
            ptr++;
        }
    }


    fill_a(mesh, a, Na, zmin, zmax);

    fastpm_smesh_add_layer_plane(mesh, xy, Nxy, a, Na);

    free(xy);
    free(a);
}

void
fastpm_smesh_add_layer_sphere(FastPMSMesh * mesh,
        double * ra, double * dec, size_t Nxy,
        double * a, size_t Na)
{
    struct FastPMSMeshLayer * layer = fastpm_smesh_add_layer_common(mesh, a, Na);

    layer->type = FASTPM_SMESH_SPHERE;
    layer->Nxy = Nxy;

    size_t i;

    layer->vec = malloc(sizeof(double) * Nxy * 3);
    layer->ra = malloc(sizeof(double) * Nxy);
    layer->dec = malloc(sizeof(double) * Nxy);

    for(i = 0; i < Nxy; i ++) {
        layer->vec[i][0] = cos(dec[i] / rad_to_degree) * cos(ra[i] / rad_to_degree);
        layer->vec[i][1] = cos(dec[i] / rad_to_degree) * sin(ra[i] / rad_to_degree);
        layer->vec[i][2] = sin(dec[i] / rad_to_degree);
        layer->ra[i] = ra[i];
        layer->dec[i] = dec[i];
    }
}

void
fastpm_smesh_add_layer_healpix(FastPMSMesh * mesh,
        int nside,
        double * a, size_t Na, MPI_Comm comm)
{
    double *ra, *dec;
    size_t npix;

    fastpm_utils_healpix_ra_dec(nside, &ra, &dec, &npix, mesh->lc->fov, comm);

    fastpm_smesh_add_layer_sphere(mesh, ra, dec, npix, a, Na);

    free(ra);
    free(dec);
}

/* automatically add healpix layers with roughly the correct
 * surface number density of mesh points */
void
fastpm_smesh_add_layers_healpix(FastPMSMesh * mesh,
        double surface_density, double volume_density,
        double amin, double amax, int maxnside,
        MPI_Comm comm)
{

    double zmin = mesh->lc->speedfactor * HorizonDistance(amax, mesh->lc->horizon);
    double zmax = mesh->lc->speedfactor * HorizonDistance(amin, mesh->lc->horizon);

    double line_density = volume_density / surface_density;

    int Na = ceil((zmax - zmin) / line_density) + 1;
    
    if(Na < 1) Na = 1; /* at least one layer */

    double * a = malloc(sizeof(double) * Na);
    size_t * nside = malloc(sizeof(size_t) * Na);

    fill_a(mesh, a, Na, zmin, zmax);

    int i;

    for(i = 0; i < Na; i ++) {
        double z = zmin + (zmax - zmin) / (Na - 1) * i;
        if(i == Na - 1) z = zmax;

        double v = sqrt(4 * M_PI / 12 * surface_density) * z;

        if(v < 1.0) v = 1.0;

        /* round to nearest power of 2 */
        nside[i] = 1L << (int64_t) (log2(v) + 0.5);

        if(maxnside > 0 && nside[i] > maxnside) {
            nside[i] = maxnside;
        }
    }

    int j = 0;
    for(i = 1; i <= Na; i ++) {
        if(nside[i] == nside[j] && i != Na) continue;
        /* nside[i] != nside[j]; j ... i - 1 (inclusive) has the same nside */
        fastpm_info("Creating smesh for Nside = %04td arange %6.4f - %6.4f\n", nside[j], a[j], a[i - 1]);
        double *ra, *dec;
        size_t npix;

        fastpm_utils_healpix_ra_dec(nside[j], &ra, &dec, &npix, mesh->lc->fov, comm);

        fastpm_smesh_add_layer_sphere(mesh, ra, dec, npix, &a[j], i - j);

        free(ra);
        free(dec);
        j = i;
    }

    free(a);
    free(nside);
}

void
fastpm_smesh_destroy(FastPMSMesh * mesh)
{
    struct FastPMSMeshLayer * layer, *next;
    for(layer = mesh->layers; layer; layer = next) {
        switch(layer->type) {
            case FASTPM_SMESH_PLANE:
                free(layer->xy);
                break;
            case FASTPM_SMESH_SPHERE:
                free(layer->ra);
                free(layer->dec);
                free(layer->vec);
            break;
        }
        free(layer->a);
        free(layer->z);
        next = layer->next;
        free(layer);
    }
    fastpm_destroy_event_handlers(&mesh->event_handlers);
    fastpm_store_destroy(mesh->last.p);
}

static void
fastpm_smesh_layer_select_active(
        FastPMSMesh * mesh,
        struct FastPMSMeshLayer * layer,
        double a0, double a1,
        FastPMStore * q
    )
{
    size_t Na = 0;

    size_t i;
    for (i = 0; i < layer->Na; i ++) {
        if(layer->a[i] >= a0 && layer->a[i] < a1) {
            Na ++;
        }
    }

    if((size_t) Na * (size_t) layer->Nxy + q->np >= q->np_upper) {
        fastpm_raise(0, "More layer points requested than np_upper (%d * %d > %td), a0 = %g a1 = %g\n",
                Na, layer->Nxy, q->np_upper, a0, a1);
    }

    size_t j = 0;
    size_t k = 0;
    size_t m=0;
    double x_temp[4];
    x_temp[3]=1;
    for(j = 0; j < layer->Nxy; j ++) {
        for(k = 0; k < layer->Na; k ++) {
            /* cast to float because aemit is float */
            float aemit = layer->a[k];
            if(aemit >= a0 && aemit < a1) {
                switch(layer->type) {
                    case FASTPM_SMESH_PLANE:
                        x_temp[0]=layer->xy[j][0];
                        x_temp[1]=layer->xy[j][1];
                        x_temp[2]=layer->z[k];
                        break;
                    case FASTPM_SMESH_SPHERE:
                        x_temp[0] = layer->vec[j][0] * layer->z[k];
                        x_temp[1] = layer->vec[j][1] * layer->z[k];
                        x_temp[2] = layer->vec[j][2] * layer->z[k];
                        break;
                }
                /* transform to simulation coordinates */
                double xo[4];
                fastpm_gldot(mesh->lc->glmatrix_inv, x_temp, xo);
                for(m = 0; m < 3; m ++) {
                    q->x[q->np][m] = xo[m];
                }
                q->aemit[q->np] = aemit;
                q->np++;
                if(q->np > q->np_upper) {
                    fastpm_raise(-1, "too many particles are created, increase np_upper!");
                }
            }
        }
    }
}

void
fastpm_smesh_select_active(FastPMSMesh * mesh,
        double a0, double a1,
        FastPMStore * q
    )
{
    struct FastPMSMeshLayer * layer;
    for(layer = mesh->layers; layer; layer = layer->next) {
        fastpm_smesh_layer_select_active(mesh, layer, a0, a1, q);
    }
}

int
fastpm_smesh_compute_potential(
        FastPMSMesh * mesh,
        PM * pm,
        FastPMGravity * gravity,
        FastPMFloat * delta_k,
        double a_f,
        double a_n
)
{
    FastPMStore p_new_now[1];
    FastPMStore p_last_now[1];

    fastpm_store_init(p_new_now, mesh->np_upper,
              PACK_POS
            | PACK_POTENTIAL
            | PACK_DENSITY
            | PACK_TIDAL
            | PACK_AEMIT,
            FASTPM_MEMORY_HEAP
    );

    fastpm_smesh_select_active(mesh, a_f, a_n, p_new_now);

    fastpm_store_wrap(p_new_now, pm->BoxSize);
    fastpm_store_decompose(p_new_now, (fastpm_store_target_func) FastPMTargetPM, pm, pm_comm(pm));

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
    
    int64_t np = p_last_now->np + p_new_now->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG, MPI_SUM, pm_comm(pm));

    if(np > 0) {
        fastpm_info("Computing potential for %ld mesh points\n", np);
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

            if(ACC[d] == PACK_DENSITY) {
                fastpm_apply_smoothing_transfer(pm, canvas, canvas, mesh->smoothing);
            }

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

    }

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

        /* transform back to observer coordinate */
        double x_temp[4], xo[4];
        int m;
        for(m = 0; m < 3; m ++) {
            x_temp[m] = p_last_then->x[i][m];
        }
        x_temp[3] = 1;

        fastpm_gldot(mesh->lc->glmatrix, x_temp, xo);

        for(m = 0; m < 3; m ++) {
            p_last_then->x[i][m] = xo[m];
        }

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

    fastpm_store_copy(p_new_now, mesh->last.p);


    fastpm_store_destroy(p_new_now);

    return 0;
}

static double
_a_to_distance(double a, void * data)
{
    void ** p = (void**) data;
    FastPMSMesh * mesh = p[0];
    double * z = p[1];
    return HorizonDistance(a, mesh->lc->horizon) * mesh->lc->speedfactor - *z;
}

static void
fill_a(FastPMSMesh * mesh, double * a, int Na,
       double zmin, double zmax)
{

    int i;
    for(i = 0; i < Na; i ++) {
        double z = zmin + (zmax - zmin) / (Na - 1) * i;
        if(i == Na - 1) z = zmax;

        void * data[2] = {mesh, &z};

        fastpm_horizon_solve(mesh->lc->horizon,
            &a[i],
            1e-7, 1,
            _a_to_distance, data);
    }

}
