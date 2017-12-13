FASTPM_BEGIN_DECLS

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMStore * p;
    /* Need a table for drift factors */

    void * gsl; // GSL solver pointer

    double glmatrix[4][4];
    double (* tileshifts)[3];
    int ntiles;
    int flatsky;

    struct {
        double * Dc;
        size_t size;
    } EventHorizonTable;

    FastPMCosmology * cosmology;

} FastPMLightCone;

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMStore * pi);

void
fastpm_lc_init(FastPMLightCone * lc,
                double speedfactor,
                double glmatrix[4][4],
                double (*tileshifts)[3], int ntiles,
                int flatsky,
                FastPMCosmology * c,
                FastPMStore * p);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
