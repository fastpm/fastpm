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

} FastPMLightCone;

typedef struct FastPMUSMesh {
    FastPMLightCone * lc;
    FastPMStore * p; /* storing the unstructred output, particles on lightcone */

    double (* tileshifts)[3];
    int ntiles;
    /* we need to apply a cut in time, because at early time we tend to write too many particles. */
    double amax; /* range for largest a; above which no particles will be written */
    double amin; /* range for smallest a; below which no particles will be written */
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
        };
        struct {
            double (* xy)[2];
        };
    };

    union {
            int Npix;
            int Nxy;
    };

    double * a;
    double * z;
    int Na;

    /* state about the last time range */
    struct {
        FastPMStore p[1];
        double a_f; /* the time potential was updated */
    } last;
    int started;

    /* Extensions */
    FastPMEventHandler * event_handlers;

} FastPMSMesh;

typedef struct FastPMLCEvent {
    FastPMEvent base;
    int is_first;
    FastPMStore * p;
    double a0;
    double a1;
} FastPMLCEvent;

void
fastpm_lc_init(FastPMLightCone * lc);

void
fastpm_usmesh_init(FastPMUSMesh * mesh,
                FastPMLightCone * lc, size_t np_upper,
                double (*tileshifts)[3], int ntiles,
                double amin, double amax);

void
fastpm_usmesh_destroy(FastPMUSMesh * mesh);

int
fastpm_usmesh_intersect(FastPMUSMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_destroy(FastPMLightCone * lc);

void
fastpm_smesh_init_plane(FastPMSMesh * mesh, FastPMLightCone * lc,
        double (*xy)[2], size_t Nxy,
        double * a, size_t Na);

void
fastpm_smesh_init_sphere(FastPMSMesh * mesh, FastPMLightCone * lc,
        double * ra, double * dec, size_t Npix,
        double * a, size_t Na);

void
fastpm_smesh_select_active(FastPMSMesh * mesh,
        double a0, double a1,
        FastPMStore * q
    );

int
fastpm_smesh_compute_potential(
        FastPMSMesh * mesh,
        PM * pm,
        FastPMGravity * gravity,
        FastPMFloat * delta_k,
        double a_f,
        double a_n);

void
fastpm_smesh_destroy(FastPMSMesh * mesh);

void fastpm_matrix_invert(double (*matrix)[4], double (*matrix_inv)[4], int matrix_size);
void healpix_ra_dec(double **ra, double **dec,long nside,long npix);

FASTPM_END_DECLS
