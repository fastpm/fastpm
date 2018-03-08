FASTPM_BEGIN_DECLS

#define FASTPM_EVENT_LC_READY "LIGHTCONEREADY"

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMCosmology * cosmology;
    FastPMHorizon * horizon;
    double speedfactor;
    double glmatrix[4][4];
    double fov; /* field of view angle. <=0 for flatsky.
                    Remember the lightcone is always along z-direction.*/

} FastPMLightCone;

typedef struct FastPMUnstructuredMesh {
    FastPMLightCone * lc;
    FastPMStore * p; /* storing the unstructred output, particles on lightcone */

    double (* tileshifts)[3];
    int ntiles;
} FastPMUnstructuredMesh;


typedef struct FastPMStructuredMesh {
    FastPMLightCone * lc;

    enum {
        FASTPM_STRUCT_MESH_SPHERE,
        FASTPM_STRUCT_MESH_PLANE,
    } type;

    union {
        struct {
            double * ra;
            double * dec;
            double (* vec)[3];
            int Npix;
        };
        struct {
            double (* xy)[2];
            int Nxy;
        };
    };

    double * z;
    int Nz;

    /* state about the last time range */
    struct {
        FastPMStore p[1];
        double z0;
        double z1;
    } last;

    /* Extensions */
    FastPMEventHandler * event_handlers;

} FastPMStructuredMesh;

typedef struct FastPMLCEvent {
    FastPMEvent base;
    FastPMStore * p;
    double z0;
    double z1;
} FastPMLCEvent;

void
fastpm_lc_init(FastPMLightCone * lc);

void
fastpm_unstruct_mesh_init(FastPMUnstructuredMesh * mesh,
                FastPMLightCone * lc, size_t np_upper,
                double (*tileshifts)[3], int ntiles);

void fastpm_unstruct_mesh_destroy(FastPMUnstructuredMesh * mesh);
int
fastpm_unstruct_mesh_intersect(FastPMUnstructuredMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_destroy(FastPMLightCone * lc);

FASTPM_END_DECLS
