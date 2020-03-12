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
    double target_volume; /* target this volume size per light cone shell during intersection search */
    /* we need to apply a cut in time, because at early time we tend to write too many particles. */
    double amax; /* range for largest a; above which no particles will be written */
    double amin; /* range for smallest a; below which no particles will be written */

    double ai; /* starting scaling factor of the current p. */
    double af; /* ending scaling factor of the current p. */
    size_t np_before; /* number of particles already written to the lightcone, on this rank.*/
    /* Extensions */
    FastPMEventHandler * event_handlers;
} FastPMUSMesh;


#define TIMESTEP_START 0
#define TIMESTEP_CUR 1
#define TIMESTEP_END 2
typedef struct FastPMLCEvent {
    FastPMEvent base;
    int whence; /* TIMESTEP_START, TIMESTEP_CUR, TIMESTEP_END for first, any, or last event. first and last shall produce no particles. */
    FastPMStore * p;
    double ai;
    double af;
} FastPMLCEvent;

void
fastpm_lc_init(FastPMLightCone * lc);

int
fastpm_lc_inside(FastPMLightCone * lc, double vec[3]);

double
fastpm_lc_distance(FastPMLightCone * lc, double x[3]);

void
fastpm_usmesh_init(FastPMUSMesh * mesh,
                FastPMLightCone * lc,
                double target_volume,
                FastPMStore * source,
                size_t np_upper,
                double (*tileshifts)[3], int ntiles,
                double amin, double amax);

void
fastpm_usmesh_destroy(FastPMUSMesh * mesh);

int
fastpm_usmesh_intersect(FastPMUSMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick,
    double a1, double a2,
    int whence, MPI_Comm comm);

void
fastpm_lc_destroy(FastPMLightCone * lc);

int
fastpm_shell_intersects_bbox(
    double xmin[3],
    double xmax[3],
    double glmatrix[4][4],
    double tileshift[3],
    double radius1,
    double radius2
);

FASTPM_END_DECLS

#endif
