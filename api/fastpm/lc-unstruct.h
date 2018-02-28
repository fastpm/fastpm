FASTPM_BEGIN_DECLS

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMCosmology * cosmology;
    FastPMHorizon * horizon;
    double speedfactor;
    double glmatrix[4][4];
    double fov; /* field of view angle. <=0 for flatsky.
                    Remember the lightcone is always along z-direction.*/

    /* private: */
    FastPMStore * unstruct; /* storing the unstructred output, particles on lightcone */

    double (* tileshifts)[3];
    int ntiles;

} FastPMLightCone;

struct FastPMStructuredMesh {
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
};

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_init(FastPMLightCone * lc, FastPMSolver * fastpm,
                double (*tileshifts)[3], int ntiles);

void
fastpm_lc_destroy(FastPMLightCone * lc);

FASTPM_END_DECLS
