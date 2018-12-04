FASTPM_BEGIN_DECLS

struct FastPMHorizon {
    FastPMCosmology * cosmology;
    size_t size;
    double da;
    double xi_a[8192];
    double growthfactor_a[8192];
    void * gsl;
};

void fastpm_horizon_init(FastPMHorizon * horizon, FastPMCosmology * cosmology);
void fastpm_horizon_destroy(FastPMHorizon * horizon);

double HorizonDistance(double a, FastPMHorizon * horizon);
double HorizonGrowthFactor(double a, FastPMHorizon * horizon);

void *
fastpm_horizon_solve_start();

void
fastpm_horizon_solve_end(void * );

int
fastpm_horizon_solve(FastPMHorizon * horizon,
    void * context,
    double * solution,
    double a_i, double a_f,
    double (*func)(double a, void * userdata),
    void * userdata);

FASTPM_END_DECLS
