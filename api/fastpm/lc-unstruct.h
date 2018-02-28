FASTPM_BEGIN_DECLS

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
        FASTPM_2DMESH_SPHERE,
        FASTPM_2DMESH_PLANE,
    } type;

    union {
        struct {
            double * ra;
            double * dec;
            int Npix;
        };
        struct {
            int Nxy;
        };
    };

    double * z;
    int Nz;
    /* private : */
    FastPMStore * q; /* storing the structred output, particles on lightcone */
} FastPMStructuredMesh;

void
fastpm_lc_init(FastPMLightCone * lc);

void
fastpm_unstruct_mesh_init(FastPMUnstructuredMesh * mesh,
                FastPMLightCone * lc, size_t np_upper,
                double (*tileshifts)[3], int ntiles);

int
fastpm_unstruct_mesh_intersect(FastPMUnstructuredMesh * mesh, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_destroy(FastPMLightCone * lc);

FASTPM_END_DECLS
