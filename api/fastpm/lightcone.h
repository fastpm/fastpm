FASTPM_BEGIN_DECLS

typedef struct {
    int compute_potential;
    /* Storage of the particles on the light cone */
    FastPMCosmology * cosmology;
    double speedfactor;
    double glmatrix[4][4];
    double fov; /* field of view angle. <=0 for flatsky.
                    Remember the lightcone is always along z-direction.*/

    /* private: */
    FastPMStore * p0; /* per step potential and tidal field.  */

    FastPMStore * p; /* storing the output, particles on lightcone */
    FastPMStore * q; /* stores the output, uniform sampling of potential on lightcone */

    double (* tileshifts)[3];
    int ntiles;
    struct {
        double * Dc;
        size_t size;
    } EventHorizonTable;

    /* Need a table for drift factors */

    void * gsl; // GSL solver pointer

} FastPMLightCone;

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMSolver * fastpm);

void
fastpm_lc_init(FastPMLightCone * lc, FastPMSolver * fastpm,
                double (*tileshifts)[3], int ntiles);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

int
fastpm_lc_compute_potential(FastPMSolver * fastpm,
        FastPMForceEvent * event,
        FastPMLightCone * lc);

FASTPM_END_DECLS
