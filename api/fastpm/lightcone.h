#ifndef FASTPM_LC_USTRUCT

#define FASTPM_LC_USTRUCT

FASTPM_BEGIN_DECLS

#define FASTPM_EVENT_LC_READY "LIGHTCONEREADY"

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMCosmology * cosmology;
    FastPMHorizon * horizon;
    double speedfactor;
    double glmatrix[4][4];
    double glmatrix_inv[4][4];
    double fov; /* field of view angle. <=0 for flatsky.
                    Remember the lightcone is always along z-direction.*/

    int octants[8]; /* 1 if the octants is to be included, enabled when fov >= 360.*/
    double tol; /* tolerance in radians in octant culling */
} FastPMLightCone;

typedef struct FastPMUSMesh {
    FastPMLightCone * lc;
    FastPMStore * source; /* source particle to monitor */
    FastPMStore * p; /* storing the unstructred output, particles on lightcone */

    double (* tileshifts)[3];
    int ntiles;
    /* we need to apply a cut in time, because at early time we tend to write too many particles. */
    double amax; /* range for largest a; above which no particles will be written */
    double amin; /* range for smallest a; below which no particles will be written */

    int is_first; /* is this the first time an event is emitted. */
    /* Extensions */
    FastPMEventHandler * event_handlers;
} FastPMUSMesh;


typedef struct FastPMLCEvent {
    FastPMEvent base;
    int is_first;
    FastPMStore * p;
    double a0;
    double a1;
} FastPMLCEvent;

void
fastpm_lc_init(FastPMLightCone * lc);

int
fastpm_lc_inside(FastPMLightCone * lc, double vec[3]);

double
fastpm_lc_distance(FastPMLightCone * lc, double x[3]);

void
fastpm_usmesh_init(FastPMUSMesh * mesh,
                FastPMLightCone * lc, FastPMStore * source,
                size_t np_upper,
                double (*tileshifts)[3], int ntiles,
                double amin, double amax);

void
fastpm_usmesh_destroy(FastPMUSMesh * mesh);

int
fastpm_usmesh_intersect(FastPMUSMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick);

void
fastpm_lc_destroy(FastPMLightCone * lc);

#include "lightcone-smesh.h"
FASTPM_END_DECLS

#endif
