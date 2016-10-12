FASTPM_BEGIN_DECLS

typedef struct {
    /* Storage of the particles on the light cone */
    FastPMStore * p;
    /* Need a table for drift factors */

    void * gsl; // GSL solver pointer

    struct {
        double * Dc;
        size_t size;
    } EventHorizonTable;

} FastPMLightCone;

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

int
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDriftFactor * drift, FastPMKickFactor * kick, FastPMStore * pi);

void
fastpm_lc_init(FastPMLightCone * lc, double speedfactor, FastPMCosmology * c, FastPMStore * p);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
