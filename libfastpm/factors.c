/*********************
 * Time intergral KDK scheme.
 * kick and drifts.
 * 
 * This code was initially modified by Jun Koda, 
 * from the original serial COLA code
 * by Svetlin Tassev.
 *
 * The kick and drift still supports a COLA compat-mode.
 * Most of the nasty factors are for COLA compat-mode
 * (not needed in PM)
 * We also added a 2LPT mode that does just 2LPT.
 *
 *  Yu Feng <rainwoodman@gmail.com> 
 *
 */

#include <math.h>
#include <string.h>
#include <assert.h>
#include <alloca.h>
#include <mpi.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "vpm.h"

static double 
Sq(double ai, double af, double aRef, double nLPT, FastPMCosmology * c, int USE_NONSTDDA);

static double 
Sphi(double ai, double af, double aRef, double nLPT, FastPMCosmology * c, int USE_NONSTDDA);

double
fastpm_solver_growth_factor(FastPMSolver * fastpm, double a)
{
    return GrowthFactor(a, fastpm->cosmology);
}

inline void
fastpm_drift_one(FastPMDriftFactor * drift, FastPMStore * p, ptrdiff_t i, double xo[3], double af)
{
    double ind;
    double dyyy, da1, da2;

    if(af == drift->af) {
        dyyy = drift->dyyy[drift->nsamples - 1];
        da1  = drift->da1[drift->nsamples - 1];
        da2  = drift->da2[drift->nsamples - 1];
    } else {
        ind = (af - drift->ai) / (drift->af - drift->ai) * (drift->nsamples - 1);
        int l = floor(ind);
        double u = l + 1 - ind;
        double v = ind - l;
        if(l + 1 >= drift->nsamples) {
            fastpm_raise(-1, "drift beyond factor's available range. ");
        }
        dyyy = drift->dyyy[l] * u + drift->dyyy[l + 1] * v;
        da1  = drift->da1[l] * u + drift->da1[l + 1] * v;
        da2  = drift->da2[l] * u + drift->da2[l + 1] * v;
    }
    int d;
    for(d = 0; d < 3; d ++) {
        double v;
        switch(drift->forcemode) {
            case FASTPM_FORCE_2LPT:
                xo[d] = p->x[i][d] + p->dx1[i][d] * da1 + p->dx2[i][d] * da2;
            break;
            case FASTPM_FORCE_ZA:
                xo[d] = p->x[i][d] + p->dx1[i][d] * da1;
            break;
            case FASTPM_FORCE_FASTPM:
            case FASTPM_FORCE_PM:
                xo[d] = p->x[i][d] + p->v[i][d] * dyyy;
            break;
            case FASTPM_FORCE_COLA:
                /* For cola, remove the lpt velocity to find the residual velocity v*/
                v = p->v[i][d] - (p->dx1[i][d]*drift->Dv1 + p->dx2[i][d]*drift->Dv2);
                xo[d] = p->x[i][d] + v * dyyy;
                xo[d] += p->dx1[i][d] * da1 + p->dx2[i][d] * da2;
            break;
        }
    }
}

inline void
fastpm_kick_one(FastPMKickFactor * kick, FastPMStore * p, ptrdiff_t i, float vo[3], double af)
{
    double ind;
    double dda, Dv1, Dv2;

    if(af == kick->af) {
        dda = kick->dda[kick->nsamples - 1];
        Dv1 = kick->Dv1[kick->nsamples - 1];
        Dv2 = kick->Dv2[kick->nsamples - 1];
    } else {
        ind = (af - kick->ai) / (kick->af - kick->ai) * (kick->nsamples - 1);
        int l = floor(ind);
        double u = l + 1 - ind;
        double v = ind - l;
        if(l + 1 >= kick->nsamples) {
            fastpm_raise(-1, "kick beyond factor's available range. ");
        }
        dda = kick->dda[l] * u + kick->dda[l + 1] * v;
        Dv1 = kick->Dv1[l] * u + kick->Dv1[l + 1] * v;
        Dv2 = kick->Dv2[l] * u + kick->Dv2[l + 1] * v;
    }

    int d;
    for(d = 0; d < 3; d++) {
        float ax = p->acc[i][d];
        if(kick->forcemode == FASTPM_FORCE_COLA) {
            ax += (p->dx1[i][d]*kick->q1 + p->dx2[i][d]*kick->q2);
        }
        vo[d] = p->v[i][d] + ax * dda;
        if(kick->forcemode == FASTPM_FORCE_COLA) {
            vo[d] += (p->dx1[i][d] * Dv1 + p->dx2[i][d] * Dv2);
        }
    }
}

// Leap frog time integration

void 
fastpm_kick_store(FastPMKickFactor * kick,
    FastPMStore * pi, FastPMStore * po, double af)
{

    if(kick->ai != pi->a_v) {
        fastpm_raise(-1, "kick is inconsitant with state.\n");
    }
    if(kick->ac != pi->a_x) {
        fastpm_raise(-1, "kick is inconsitant with state.\n");
    }
    int np = pi->np;

    // Kick using acceleration at a= ac
    // Assume forces at a=ac is in particles->force

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        float vo[3];
        fastpm_kick_one(kick, pi, i, vo, af);
        for(d = 0; d < 3; d++) {
            po->v[i][d] = vo[d];
        }
    }

    //velocity is now at a= avel1
    po->a_v = af;
}

static double G_p(double a, FastPMCosmology * c)
{
    /* integral of G_p */
    return GrowthFactor(a, c);
}
static double g_p(double a, FastPMCosmology * c)
{
    return DGrowthFactorDa(a, c);
}

static double G_f(double a, FastPMCosmology * c)
{
    /* integral of g_f */
    return a * a * a * HubbleEa(a, c) * g_p(a, c);
}

static double g_f(double a, FastPMCosmology * c)
{
    double dDda = DGrowthFactorDa(a, c);
    double E = HubbleEa(a, c);
    double dEda = DHubbleEaDa(a, c);
    double d2Dda2 = D2GrowthFactorDa2(a, c);

    double g_f = 3 * a * a * E * dDda
                   + a * a * a * dEda * dDda
                   + a * a * a * E * d2Dda2;
    return g_f;
}
void fastpm_kick_init(FastPMKickFactor * kick, FastPMSolver * fastpm, double ai, double ac, double af)
{
    FastPMCosmology * c = fastpm->cosmology;
    kick->forcemode = fastpm->config->FORCE_TYPE;

    double OmegaM = c->OmegaM;
    double Om143 = pow(OmegaA(ac, c), 1.0/143.0);
    double growth1 = GrowthFactor(ac, c);

    kick->q1 = growth1;
    kick->q2 = growth1*growth1*(1.0 + 7.0/3.0*Om143);

    kick->nsamples = 32;
    int i;

    double Dv1i = GrowthFactor(ai, c) * ai * ai * HubbleEa(ai, c) * DLogGrowthFactor(ai, c);
    double Dv2i = GrowthFactor2(ai, c) * ai * ai * HubbleEa(ai, c) * DLogGrowthFactor2(ai, c);
    for(i = 0; i < kick->nsamples; i ++) {
        double ae = ai * (1.0 * (kick->nsamples - 1 - i) / (kick->nsamples - 1))
                  + af * (1.0 * i / (kick->nsamples - 1));

        if(kick->forcemode == FASTPM_FORCE_FASTPM) {
            kick->dda[i] = -1.5 * OmegaM
               * 1 / (ac * ac * HubbleEa(ac, c))
               * (G_f(ae, c) - G_f(ai, c)) / g_f(ac, c);
        } else {
            kick->dda[i] = -1.5 * OmegaM * Sphi(ai, ae, ac, fastpm->config->nLPT, c, kick->forcemode == FASTPM_FORCE_COLA);
        }
        kick->Dv1[i] = GrowthFactor(ae, c) * ae * ae * HubbleEa(ae, c) * DLogGrowthFactor(ae, c) - Dv1i;
        kick->Dv2[i] = GrowthFactor2(ae, c) * ae * ae * HubbleEa(ae, c) * DLogGrowthFactor2(ae, c) - Dv2i;
    }

    kick->ai = ai;
    kick->ac = ac;
    kick->af = af;

}

void
fastpm_drift_init(FastPMDriftFactor * drift, FastPMSolver * fastpm,
                double ai, double ac, double af)
{
    FastPMCosmology * c = fastpm->cosmology;
    drift->forcemode = fastpm->config->FORCE_TYPE;

    drift->nsamples = 32;
    int i;

    for(i = 0; i < drift->nsamples; i ++ ) {
        double ae = ai * (1.0 * (drift->nsamples - 1 - i) / (drift->nsamples - 1))
                  + af * (1.0 * i / (drift->nsamples - 1));

        if(drift->forcemode == FASTPM_FORCE_FASTPM) {
            drift->dyyy[i] = 1 / (ac * ac * ac * HubbleEa(ac, c))
                        * (G_p(ae, c) - G_p(ai, c)) / g_p(ac, c);
        } else {
            drift->dyyy[i] = Sq(ai, ae, ac, fastpm->config->nLPT, c, drift->forcemode == FASTPM_FORCE_COLA);
        }
        drift->da1[i] = GrowthFactor(ae, c) - GrowthFactor(ai, c);    // change in D_1lpt
        drift->da2[i] = GrowthFactor2(ae, c) - GrowthFactor2(ai, c);  // change in D_2lpt
    }
    drift->af = af;
    drift->ai = ai;
    drift->ac = ac;
    drift->Dv1 = GrowthFactor(ac, c) * ac * ac * HubbleEa(ac, c) * DLogGrowthFactor(ac, c);
    drift->Dv2 = GrowthFactor2(ac, c) * ac * ac * HubbleEa(ac, c) * DLogGrowthFactor2(ac, c);
}

void
fastpm_drift_store(FastPMDriftFactor * drift,
               FastPMStore * pi, FastPMStore * po,
               double af)
{
    if(drift->ai != pi->a_x) {
        fastpm_raise(-1, "drift is inconsitant with state.\n");
    }
    if(drift->ac != pi->a_v) {
        fastpm_raise(-1, "drift is inconsitant with state.\n");
    }
    int np = pi->np;

    int i;
    // Drift
#pragma omp parallel for
    for(i=0; i<np; i++) {
        double xo[3] = {0};
        fastpm_drift_one(drift, pi, i, xo, drift->af);
        int d;
        for(d = 0; d < 3; d ++) {
            po->x[i][d] = xo[d];
        }
    }
    po->a_x = af;
}

//
// Functions for our modified time-stepping (used when StdDA=0):
//
struct iparam {
    FastPMCosmology * cosmology;
    double nLPT;
};

double gpQ(double a, double nLPT) { 
    return pow(a, nLPT);
}

static double stddriftfunc (double a, struct iparam * iparam) {
    return 1 / (pow(a, 3) * HubbleEa(a, iparam->cosmology));
}

static double nonstddriftfunc (double a, struct iparam * iparam) {
    return gpQ(a, iparam->nLPT)/(pow(a, 3) * HubbleEa(a, iparam->cosmology));
}

static double stdkickfunc (double a, struct iparam * iparam) {
    return 1/ (pow(a, 2) * HubbleEa(a, iparam->cosmology));
}

static double integrand(double a, void * params) {
    void ** p = (void**) params;
    double (*func)(double a, struct iparam * s) = p[0];
    struct iparam * s = p[1];
    return func(a, s);
}

double integrate(double ai, double af,
        struct iparam * iparam,
        double (*func)(double , struct iparam * )) {

    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (5000);

    gsl_function F;
    double error;
    double result;

    F.params = (void*[]){func, iparam};
    F.function = integrand;

    gsl_integration_qag (&F, ai, af, 0, 1e-8, 5000, 6,
            w, &result, &error); 

    gsl_integration_workspace_free (w);
    return result;
}

/*     
       When StdDA=0, one needs to set nLPT.
       assumes time dep. for velocity = B a^nLPT
       nLPT is a real number. Sane values lie in the range (-4,3.5). Cannot be 0, but of course can be -> 0 (say 0.001).
       See Section A.3 of TZE.
       */

static double 
Sq(double ai, double af, double aRef, double nLPT, FastPMCosmology * c, int USE_NONSTDDA)
{
    double resultstd, result;
    struct iparam iparam[1];
    iparam->cosmology = c;
    iparam->nLPT = nLPT;

    resultstd = integrate(ai, af, iparam, stddriftfunc);

    result = integrate(ai, af, iparam, nonstddriftfunc);
    result /= gpQ(aRef, nLPT);

    /*
    fastpm_info("ref time = %6.4f, std drift =%g, non std drift = %g \n",
        aRef, resultstd, result); */

    if (USE_NONSTDDA)
        return result;
    else
        return resultstd;
}

double DERgpQ(double a, double nLPT) { 
    /* This must return d(gpQ)/da */
    return nLPT*pow(a, nLPT-1);
}


static double 
Sphi(double ai, double af, double aRef, double nLPT, FastPMCosmology * c, int USE_NONSTDDA)
{
    double result;
    double resultstd;

    struct iparam iparam[1];
    iparam->cosmology = c;
    iparam->nLPT = nLPT;

    result = (gpQ(af, nLPT) - gpQ(ai, nLPT)) * aRef 
        / (pow(aRef, 3) * HubbleEa(aRef, c) * DERgpQ(aRef, nLPT));

    resultstd = integrate(ai, af, iparam, stdkickfunc);

    /*
    fastpm_info("ref time = %6.4f, std kick = %g, non std kick = %g\n",
            aRef, resultstd, result); */

    if (USE_NONSTDDA) {
        return result;
    } else {
        return resultstd;
    }
}
