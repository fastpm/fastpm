#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#define STEF_BOLT 2.85087e-48  // in units: [ h * (10^10Msun/h) * s^-3 * K^-4 ]
#define rho_crit 27.7455       // rho_crit0 in mass/length
#define LIGHT 9.715614e-15     // in units: [ h * (Mpc/h) * s^-1 ]

#define kB 8.617333262145e-5   // boltzman in eV/K

double HubbleDistance = 2997.92458; /* Mpc/h */
double HubbleConstant = 100.0;      /* km/s / Mpc/h */


void
fastpm_cosmology_init(FastPMCosmology * c)
{
    /* prepare the interpolation object for FD tables. */
    if (c->N_ncdm > 0) {
        FastPMFDInterp * FDinterp = malloc(sizeof(FDinterp[0]));
        fastpm_fd_interp_init(FDinterp);
        c->FDinterp = FDinterp;
    }

    // Compute Omega_cdm assuming all ncdm is matter like
    c->Omega_cdm = c->Omega_m - Omega_ncdmTimesHubbleEaSq(1, c);

    // Set Omega_Lambda at z=0 to give no curvature
    c->Omega_Lambda = 1 - c->Omega_m - Omega_r(c);

}

void
fastpm_cosmology_destroy(FastPMCosmology * c)
{
    if (c->N_ncdm > 0) {
        fastpm_fd_interp_destroy(c->FDinterp);
        free(c->FDinterp);
    }
}

double Omega_g(FastPMCosmology * c)
{
    /* FIXME: This is really Omega_g0. Might want to redefine all Omegas in the code */
    return 4 * STEF_BOLT * pow(c->T_cmb, 4) / pow(LIGHT, 3) / rho_crit / pow(c->h, 2);
}

double Gamma_nu(FastPMCosmology * c)
{
    /* nu to photon temp ratio today */
    if (c->N_nu == 0) {
        return 0;
    } else {
        return pow(4./11., 1./3.) * pow(c->N_eff / c->N_nu, 1./4.);
    }
}

double Omega_ur(FastPMCosmology * c)
{
    /* Omega_ur0. This is the energy density of all massless nus */
    int N_ur = c->N_nu - c->N_ncdm;    // number of massless nus (different to CLASS defn)
    return 7./8. * N_ur * pow(Gamma_nu(c), 4) * Omega_g(c);
}

double Omega_r(FastPMCosmology * c)
{
    /* Omega_r0. This is the energy density of all radiation-like particles.
       FIXME: really this is Omega_g+ur, not r, because it doesn't have ncdm,r in it 
       and we define m with ncdm,m. this doesn't affect rest of code, but maybe redefine.*/
    return Omega_g(c) + Omega_ur(c);
}

static double getFtable(int F_id, double y, FastPMCosmology * c)
{
    /* Not using size as an arg, it's globally defined
       Gets the interpolated value of Ftable[F_id] at y
       F_id: 1 for F, 2 for F', 3 for F'' */
    return fastpm_do_fd_interp(c->FDinterp, F_id, y);
}

static double Fconst(int ncdm_id, FastPMCosmology * c)
{
    /* This is a cosmology dependent constant 
       which is the argument divided by a of F, DF, DDF */
    double T_nu = Gamma_nu(c) * c->T_cmb;
    return c->m_ncdm[ncdm_id] / (kB * T_nu);
}

double Omega_ncdm_iTimesHubbleEaSq(double a, int ncdm_id, FastPMCosmology * c)   
{
    /* Omega_ncdm_i(a) * E(a)^2 */
    
    double A = 15. / pow(M_PI, 4) * pow(Gamma_nu(c), 4) * Omega_g(c);
    double Fc = Fconst(ncdm_id, c);
    double F = getFtable(1, Fc*a, c);         //row 1 for F
    
    return A / (a*a*a*a) * F;
}

double Omega_ncdmTimesHubbleEaSq(double a, FastPMCosmology * c)   
{
    // sum over ncdm species
    double res = 0;
    for (int i=0; i<c->N_ncdm; i++) {
        res += Omega_ncdm_iTimesHubbleEaSq(a, i, c);
    }
    return res;
}

double DOmega_ncdmTimesHubbleEaSqDa(double a, FastPMCosmology * c)
{
    double A = 15. / pow(M_PI, 4) * pow(Gamma_nu(c), 4) * Omega_g(c);
    
    double OncdmESq = Omega_ncdmTimesHubbleEaSq(a,c);
    
    double FcDF = 0;
    for (int i=0; i<c->N_ncdm; i++) {
        double Fc = Fconst(i, c);
        double DF = getFtable(2, Fc*a, c);
        FcDF += Fc * DF;    //row 2 for F'
    }
    
    return -4. / a * OncdmESq + A / (a*a*a*a) * FcDF;
}

double D2Omega_ncdmTimesHubbleEaSqDa2(double a, FastPMCosmology * c)
{
    double A = 15. / pow(M_PI, 4) * pow(Gamma_nu(c), 4) * Omega_g(c);
    
    double OncdmESq = Omega_ncdmTimesHubbleEaSq(a,c);
    double DOncdmESqDa = DOmega_ncdmTimesHubbleEaSqDa(a,c);
    
    double FcFcDDF = 0;
    for (int i=0; i<c->N_ncdm; i++) {
        double Fc = Fconst(i, c);
        double DDF = getFtable(3, Fc*a, c);
        FcFcDDF += Fc * Fc * DDF;    //row 3 for F''
    }
    
    return -12. / (a*a) * OncdmESq - 8. / a * DOncdmESqDa + A / (a*a*a*a) * FcFcDDF;
}

double HubbleEa(double a, FastPMCosmology * c)
{
    /* H(a) / H0 
       ncdm is NOT assumed to be matter like here */
    return sqrt(Omega_r(c) / (a*a*a*a) + c->Omega_cdm / (a*a*a) + Omega_ncdmTimesHubbleEaSq(a, c) + c->Omega_Lambda);
}

double Omega_cdm_a(double a, FastPMCosmology * c)
{
    /* as a func of a 
       FIXME: remove the _a and put in 0s for today's value elsewhere */
    double E = HubbleEa(a, c);
    return c->Omega_cdm / (a*a*a) / (E*E);
}

double Omega_m(double a, FastPMCosmology * c){
    /* Total matter component (cdm + ncdm) assming all ncdm is matter like */
    double E = HubbleEa(a, c);
    return c->Omega_m / (a*a*a) / (E*E);
}

double DHubbleEaDa(double a, FastPMCosmology * c)
{
    /* d E / d a*/
    double E = HubbleEa(a, c);
    double DOncdmESqDa = DOmega_ncdmTimesHubbleEaSqDa(a,c);
    
    return 0.5 / E * ( - 4 * Omega_r(c) / pow(a,5) - 3 * c->Omega_cdm / pow(a,4) + DOncdmESqDa );
}

double D2HubbleEaDa2(double a, FastPMCosmology * c)
{
    double E = HubbleEa(a,c);
    double dEda = DHubbleEaDa(a,c);
    double D2OncdmESqDa2 = D2Omega_ncdmTimesHubbleEaSqDa2(a,c);

    return 0.5 / E * ( 20 * Omega_r(c) / pow(a,6) + 12 * c->Omega_cdm / pow(a,5) + D2OncdmESqDa2 - 2 * pow(dEda,2) );
}

static double growth_int(double a, void *param)
{
    double * p = (double*) param;
    double Omega_m = p[0];
    double Omega_Lambda = p[1];
    return pow(a / (Omega_m + (1 - Omega_m - Omega_Lambda) * a + Omega_Lambda * a * a * a), 1.5);
}


static double solve_growth_int(double a, FastPMCosmology * c)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {c->Omega_m, c->Omega_Lambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41,
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return HubbleEa(a, c) * result;
}

static int growth_ode(double a, const double y[], double dyda[], void *params)
{
    FastPMCosmology* c = (FastPMCosmology * ) params;
    
    const double E = HubbleEa(a, c);
    const double dEda = DHubbleEaDa(a, c);
    
    double dydlna[4];
    dydlna[0] = y[1];
    dydlna[1] = - (2. + a / E * dEda) * y[1] + 1.5 * Omega_m(a, c) * y[0];
    dydlna[2] = y[3];
    dydlna[3] = - (2. + a / E * dEda) * y[3] + 1.5 * Omega_m(a, c) * (y[2] - y[0]*y[0]);
    
    //divide by  a to get dyda
    for (int i=0; i<4; i++){
        dyda[i] = dydlna[i] / a;
    }
    
    return GSL_SUCCESS;
}

static ode_soln growth_ode_solve(double a, FastPMCosmology * c)
{
    /* This returns an array of {d1, F1, d2, F2} (unnormalised) */
    gsl_odeiv2_system F;
    F.function = &growth_ode;
    F.jacobian = NULL;
    F.dimension = 4;
    F.params = (void*) c;
    
    gsl_odeiv2_driver * drive 
        = gsl_odeiv2_driver_alloc_standard_new(&F,
                                               gsl_odeiv2_step_rkf45, 
                                               1e-6,
                                               1e-8,
                                               1e-8,
                                               1,
                                               1);
    
    // assume matter domination
    double aini = 0.00625;  // FIXME: need to make sure this is less than the starting a. For now using z=159.
    double yini[4];
    yini[0] = aini;
    yini[1] = aini;
    yini[2] = - 3./7. * aini*aini;
    yini[3] = 2 * yini[2];
    
    int status = gsl_odeiv2_driver_apply(drive, &aini, a, yini);
    if (status != GSL_SUCCESS) {
        fastpm_raise(-1, "Growth ODE unsuccesful at a=%g.", a);
    }
    
    gsl_odeiv2_driver_free(drive);
    
    ode_soln soln;
    soln.y0 = yini[0];
    soln.y1 = yini[1];
    soln.y2 = yini[2];
    soln.y3 = yini[3];
    
    return soln;
}

void fastpm_growth_info_init(FastPMGrowthInfo * growth_info, double a, FastPMCosmology * c) {
    growth_info->a = a;
    growth_info->c = c;

    switch (c->growth_mode) {
        case FASTPM_GROWTH_MODE_LCDM: {
            double d1 = solve_growth_int(a, c);
            double d1_a1 = solve_growth_int(1, c);
            double Om = Omega_m(a, c);

            growth_info->D1 = d1 / d1_a1;
            growth_info->f1 = pow(Om, 5./9.);
            growth_info->D2 = growth_info->D1 * growth_info->D1 * pow(Om/Omega_m(1, c), -1./143.);
            growth_info->f2 = 2 * pow(Om, 6./11.);
        break; }
        case FASTPM_GROWTH_MODE_ODE: {
            ode_soln soln = growth_ode_solve(a, c);
            ode_soln soln_a1 = growth_ode_solve(1, c);

            growth_info->D1 = soln.y0 / soln_a1.y0;
            growth_info->f1 = soln.y1 / soln.y0;    /* f = d log D / d log a. Note soln.y1 is d d1 / d log a */
            growth_info->D2 = soln.y2 / soln_a1.y2;
            growth_info->f2 = soln.y3 / soln.y2;
        break; }
        default:
            fastpm_raise(-1, "Please enter a valid growth mode.\n");
    }
}

double DGrowthFactorDa(FastPMGrowthInfo * growth_info) {
    /* dD/da */
    /* FIXME: Technically the ODE version should agree with LCDM version
              for Lambda+CDM background,
              but to ensure backwards compatibility both versions are here. */
    double a = growth_info->a;
    FastPMCosmology * c = growth_info->c;
    
    double ans = 0;
    switch (c->growth_mode) {
        case FASTPM_GROWTH_MODE_LCDM: {
            double E = HubbleEa(a, c);
            double EI = solve_growth_int(1.0, c);
            double t1 = DHubbleEaDa(a, c) * growth_info->D1 / E;
            double t2 = E * pow(a * E, -3) / EI;
            ans = t1 + t2;
        break; }
        case FASTPM_GROWTH_MODE_ODE:
            ans = growth_info->f1 * growth_info->D1 / growth_info->a;
        break;
        default:
            fastpm_raise(-1, "Please enter a valid growth mode.\n");
    }
    return ans;
}

double D2GrowthFactorDa2(FastPMGrowthInfo * growth_info) {
    /* d^2 D1 / da^2 */
    /* FIXME: Technically the ODE version should agree with LCDM version
              for Lambda+CDM background,
              but to ensure backwards compatibility both versions are here. */
    double a = growth_info->a;
    FastPMCosmology * c = growth_info->c;
    
    double ans = 0;
    switch (c->growth_mode) {
        case FASTPM_GROWTH_MODE_LCDM: {
            double d2Eda2 = D2HubbleEaDa2(a, c);
            double dEda = DHubbleEaDa(a, c);
            double E = HubbleEa(a, c);
            double EI = solve_growth_int(1., c);
            double t1 = d2Eda2 * growth_info->D1 / E;
            double t2 = (dEda + 3 / a * E) * pow(a * E, -3) / EI;
            ans = t1 - t2;
        break; }
        case FASTPM_GROWTH_MODE_ODE: {
            double E = HubbleEa(a, c);
            double dEda = DHubbleEaDa(a, c);
            double D1 = growth_info->D1;
            double f1 = growth_info->f1;

            ans -= (3. + a / E * dEda) * f1;
            ans += 1.5 * Omega_m(a, c);
            ans *= D1 / (a*a);
        break; }
        default:
            fastpm_raise(-1, "Please enter a valid growth mode.\n");
    }
    return ans;
}

static double
comoving_distance_int(double a, void * params)
{
    FastPMCosmology * c = (FastPMCosmology * ) params;
    return 1. / (a * a * HubbleEa(a, c));
}

double ComovingDistance(double a, FastPMCosmology * c) {

    /* We tested using ln_a doesn't seem to improve accuracy */
    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &comoving_distance_int;
    F.params = (void*) c;

    gsl_integration_qag(&F, a, 1., 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41,
            workspace, &result, &abserr); 
            //lowered tol by /10 to avoid error from round off (maybe I need to make my paras more accurate)

    gsl_integration_workspace_free(workspace);

    return result;
}

#ifdef TEST_COSMOLOGY
int main() {
    /* the old COLA growthDtemp is 6 * pow(1 - c.OmegaM, 1.5) times growth */
    double a;
    FastPMCosmology c[1] = {{
        .OmegaM = 0.3,
        .OmegaLambda = 0.7
    }};

    printf("OmegaM D dD/da d2D/da2 D2 E dE/dA d2E/da2 \n");
    for(c->OmegaM = 0.1; c->OmegaM < 0.6; c->OmegaM += 0.1) {
        double f = 6 * pow(1 - c->OmegaM, 1.5);
        c->OmegaLambda = 1 - c->OmegaM;
        double a = 0.8;
        printf("%g %g %g %g %g %g %g %g %g\n",
            c->OmegaM, 
            ComovingDistance(a, c),
            GrowthFactor(a, c),
            DGrowthFactorDa(a, c),
            D2GrowthFactorDa2(a, c),

            GrowthFactor2(a, c),
            HubbleEa(a, c),
            DHubbleEaDa(a, c),
            D2HubbleEaDa2(a, c)
            );
    }
}


#endif
