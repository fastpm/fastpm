FASTPM_BEGIN_DECLS

extern double HubbleConstant;
extern double HubbleDistance;

typedef enum {
    FASTPM_GROWTH_MODE_LCDM = 0,
    FASTPM_GROWTH_MODE_ODE = 1,
} FastPMGrowthMode;

struct FastPMCosmology {
    double h;
    double Omega_m;
    double Omega_cdm;
    double Omega_ncdm;
    double Omega_k;
    double Omega_Lambda;     // Omega of dark energy at z=0
    double w0;
    double wa;
    double T_cmb;
    //double T_nu;           // todays neutrino temperature HARD CODED FOR NOW
    double N_eff;
    int N_nu;                // total number of neutrino species (massive and massless)
    double m_ncdm[3];        // masses of massive neutrinos (ncdm) for now assume max of 3 ncdms
    int N_ncdm;
    int ncdm_freestreaming;  // bool: treat ncdm as free-streaming?
    int ncdm_matterlike;     // bool: treat ncdm as matter-like?

    FastPMGrowthMode growth_mode;
    FastPMFDInterp * FDinterp;
};

double interpolate(const double xa[], const double ya[], size_t size, double xi);

double Omega_g(FastPMCosmology * c);
double Gamma_nu(FastPMCosmology * c);
double Omega_ur(FastPMCosmology * c);
double Omega_r(FastPMCosmology * c);

double Omega_ncdm_iTimesHubbleEaSq(double a, int ncdm_id, FastPMCosmology * c);

double Omega_ncdmTimesHubbleEaSq(double a, FastPMCosmology * c);
double DOmega_ncdmTimesHubbleEaSqDa(double a, FastPMCosmology * c);
double D2Omega_ncdmTimesHubbleEaSqDa2(double a, FastPMCosmology * c);

double Omega_DE_TimesHubbleEaSq(double a, FastPMCosmology * c);
double DOmega_DE_TimesHubbleEaSqDa(double a, FastPMCosmology * c);
double D2Omega_DE_TimesHubbleEaSqDa2(double a, FastPMCosmology * c);

double HubbleEa(double a, FastPMCosmology * c);
double Omega_cdm_a(double a, FastPMCosmology * c);
double Omega_m(double a, FastPMCosmology * c);
double Omega_source(double a, FastPMCosmology * c);
double DHubbleEaDa(double a, FastPMCosmology * c);
double D2HubbleEaDa2(double a, FastPMCosmology * c);

typedef struct {
    double y0;
    double y1;
    double y2;
    double y3;
} ode_soln;

typedef struct FastPMGrowthInfo {
    /* Object to store solutions to growth ode at
       a certain scale factor, a, and for a certain
       cosmology, c. */
    double a;
    FastPMCosmology * c;
    double D1;   /* growth factor normalised to 1 today */
    double D2;
    double f1;   /* dlogD1 / dloga */
    double f2;
} FastPMGrowthInfo;

void fastpm_cosmology_init(FastPMCosmology * c);
void fastpm_cosmology_destroy(FastPMCosmology * c);

void fastpm_growth_info_init(FastPMGrowthInfo * growth_info, double a, FastPMCosmology * c);

double DGrowthFactorDa(FastPMGrowthInfo * growth_info);
double D2GrowthFactorDa2(FastPMGrowthInfo * growth_info);

double ComovingDistance(double a, FastPMCosmology * c);
double OmegaA(double a, FastPMCosmology * c);

FASTPM_END_DECLS
