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

typedef struct FastPMUSMesh {
    FastPMLightCone * lc;
    FastPMStore * p; /* storing the unstructred output, particles on lightcone */

    double (* tileshifts)[3];
    int ntiles;
} FastPMUSMesh;


typedef struct FastPMSMesh {
    FastPMLightCone * lc;

    enum {
        FASTPM_SMESH_SPHERE,
        FASTPM_SMESH_PLANE,
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

} FastPMSMesh;

typedef struct FastPMLCEvent {
    FastPMEvent base;
    FastPMStore * p;
    double z0;
    double z1;
} FastPMLCEvent;

void
fastpm_lc_init(FastPMLightCone * lc);

void
fastpm_usmesh_init(FastPMUSMesh * mesh,
                FastPMLightCone * lc, size_t np_upper,
                double (*tileshifts)[3], int ntiles);

void fastpm_usmesh_destroy(FastPMUSMesh * mesh);
int
fastpm_usmesh_intersect(FastPMUSMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_destroy(FastPMLightCone * lc);

FASTPM_END_DECLS
