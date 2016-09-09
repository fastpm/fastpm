FASTPM_BEGIN_DECLS

typedef struct {
    /* remember the solver, might be useful. */
    FastPMSolver * fastpm;
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
fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMKick * kick, FastPMStore * pi);

void
fastpm_lc_init(FastPMLightCone * lc, double speedfactor, FastPMSolver * fastpm, size_t np_upper);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
