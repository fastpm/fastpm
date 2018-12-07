FASTPM_BEGIN_DECLS

extern double HubbleConstant;
extern double HubbleDistance;

struct FastPMCosmology {
    double h;
    double Omega_cdm;
    double Omega_Lambda;
    double T_cmb;    /*related to omegaR*/
    double N_eff;  //N_ur;      /*this is N_eff*/ //actually i might just make this number o fmassless neutrinos, just adds to rad
    double M_nu[3]; ///    m_ncdm[3]; for now assume 3 nus of same mass
    int N_nu;  //N_ncdm; must be less than 3 + 1.
};

double interpolate(const double xa[], const double ya[], size_t size, double xi);

double Omega_g(FastPMCosmology * c);

double getFtable(int F_id, double y);
double Fconst(FastPMCosmology * c);

double OmegaNuTimesHubbleEaSq(double a, FastPMCosmology * c);
double DOmegaNuTimesHubbleEaSqDa(double a, FastPMCosmology * c);
double D2OmegaNuTimesHubbleEaSqDa2(double a, FastPMCosmology * c);

double HubbleEa(double a, FastPMCosmology * c);
double Omega_cdm_a(double a, FastPMCosmology * c);
double OmegaA(double a, FastPMCosmology * c);
double DHubbleEaDa(double a, FastPMCosmology * c);
double D2HubbleEaDa2(double a, FastPMCosmology * c);

double OmegaSum(double a, FastPMCosmology* c);

typedef struct {
    double y0;
    double y1;
    double y2;
    double y3;
} ode_soln;

//dont put these in .h
//static int growth_odeNu(double a, const double y[], double dyda[], void *params);
//static ode_soln growth_ode_solve(double a, FastPMCosmology * c);

double growth(double a, FastPMCosmology * c);
double DgrowthDlna(double a, FastPMCosmology * c);
double growth2(double a, FastPMCosmology * c);
double Dgrowth2Dlna(double a, FastPMCosmology * c);

double GrowthFactor(double a, FastPMCosmology * c);
double GrowthFactor2(double a, FastPMCosmology * c);

double DLogGrowthFactor(double a, FastPMCosmology * c);
double DLogGrowthFactor2(double a, FastPMCosmology * c);

double DGrowthFactorDa(double a, FastPMCosmology * c);
double D2GrowthFactorDa2(double a, FastPMCosmology * c);

double ComovingDistance(double a, FastPMCosmology * c);
double OmegaA(double a, FastPMCosmology * c);

FASTPM_END_DECLS
