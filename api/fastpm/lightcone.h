FASTPM_BEGIN_DECLS

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMCosmology * cosmology;
    double speedfactor;
    double glmatrix[4][4];
    int flatsky;

    /* private: */
    FastPMStore * p;
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
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMStore * pi);

void
fastpm_lc_init(FastPMLightCone * lc, FastPMStore * p,
                double (*tileshifts)[3], int ntiles);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
