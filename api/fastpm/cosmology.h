FASTPM_BEGIN_DECLS

extern double HubbleConstant;
extern double HubbleDistance;

struct FastPMCosmology {
    double OmegaM;
    double OmegaLambda;
};

double GrowthFactor(double a, FastPMCosmology * c);
double GrowthFactor2(double a, FastPMCosmology * c);

double DLogGrowthFactor(double a, FastPMCosmology * c);
double DLogGrowthFactor2(double a, FastPMCosmology * c);

double HubbleEa(double a, FastPMCosmology * c);
double DHubbleEaDa(double a, FastPMCosmology * c);
double D2HubbleEaDa2(double a, FastPMCosmology * c);
double DGrowthFactorDa(double a, FastPMCosmology * c);
double D2GrowthFactorDa2(double a, FastPMCosmology * c);

double ComovingDistance(double a, FastPMCosmology * c);
double OmegaA(double a, FastPMCosmology * c);


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
