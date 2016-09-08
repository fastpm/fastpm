FASTPM_BEGIN_DECLS

typedef struct {
    double LightSpeedFactor;

    /* Storage of the particles on the light cone */
    FastPMStore * p;
    /* Need a table for drift factors */

	void *s; // GSL solver pointer

    struct {
        double * Dc;
        size_t size;
    } EventHorizonTable;

} FastPMLightCone;

struct funct_params { FastPMLightCone *lc; double a, b;};

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

int
fastpm_lc_intersect(FastPMLightCone * lc, double * solution, double a, double b);

void
fastpm_lc_init(FastPMLightCone * lc, Cosmology CP, size_t np_upper);

void
fastpm_lc_destroy(FastPMLightCone * lc);

double
fastpm_lc_horizon(FastPMLightCone * lc, double a);

FASTPM_END_DECLS
