FASTPM_BEGIN_DECLS

typedef struct {
    double LightSpeedFactor;

    /* storage of the particles on light cone */
    FastPMStore * p;
    /* need a table for drift factors */

    struct {
        double * Dc;
        size_t size;
    } EventHorizonTable;

} FastPMLightCone;

void
fastpm_lc_init(FastPMLightCone * lc, Cosmology CP, size_t np_upper);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
