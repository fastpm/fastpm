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

static inline void
fastpm_drift_lookup(FastPMDriftFactor * drift, double af, double * dyyy, double * da1, double * da2)
{
    double ind;

    if(af == drift->af) {
        *dyyy = drift->dyyy[drift->nsamples - 1];
        *da1  = drift->da1[drift->nsamples - 1];
        *da2  = drift->da2[drift->nsamples - 1];
        return;
    }
    if(af == drift->ai) {
        *dyyy = drift->dyyy[0];
        *da1  = drift->da1[0];
        *da2  = drift->da2[0];
        return;
    }
    {
        ind = (af - drift->ai) / (drift->af - drift->ai) * (drift->nsamples - 1);
        int l = floor(ind);
        double u = l + 1 - ind;
        double v = ind - l;
        if(l + 1 >= drift->nsamples) {
            fastpm_raise(-1, "drift beyond factor's available range. ");
        }
        *dyyy = drift->dyyy[l] * u + drift->dyyy[l + 1] * v;
        *da1  = drift->da1[l] * u + drift->da1[l + 1] * v;
        *da2  = drift->da2[l] * u + drift->da2[l + 1] * v;
    }
}

inline void
fastpm_drift_one(FastPMDriftFactor * drift, FastPMStore * p, ptrdiff_t i, double xo[3], double af)
{

    double dyyy_f, da1_f, da2_f;
    double dyyy_i, da1_i, da2_i;
    double dyyy, da1, da2;

    fastpm_drift_lookup(drift, af, &dyyy_f, &da1_f, &da2_f);
    fastpm_drift_lookup(drift, p->meta.a_x, &dyyy_i, &da1_i, &da2_i);

    dyyy = dyyy_f - dyyy_i;
    da1 = da1_f - da1_i;
    da2 = da2_f - da2_i;

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
        /* if PGDCorrection is enabled, add it */
        if(p->pgdc) {
            /* no drift; to protect the pgdc line */
            if (drift->ai == drift->af) continue;
            xo[d] += 0.5 * p->pgdc[i][d] * dyyy / drift->dyyy[drift->nsamples-1];
        }
    }
}
static inline void
fastpm_kick_lookup(FastPMKickFactor * kick, double af, double * dda, double * Dv1, double * Dv2)
{
    double ind;

    if(af == kick->af) {
        *dda = kick->dda[kick->nsamples - 1];
        *Dv1 = kick->Dv1[kick->nsamples - 1];
        *Dv2 = kick->Dv2[kick->nsamples - 1];
        return;
    }
    if(af == kick->ai) {
        *dda = kick->dda[0];
        *Dv1 = kick->Dv1[0];
        *Dv2 = kick->Dv2[0];
        return;
    }
    {
        ind = (af - kick->ai) / (kick->af - kick->ai) * (kick->nsamples - 1);
        int l = floor(ind);
        double u = l + 1 - ind;
        double v = ind - l;
        if(l + 1 >= kick->nsamples) {
            fastpm_raise(-1, "kick beyond factor's available range. ");
        }
        *dda = kick->dda[l] * u + kick->dda[l + 1] * v;
        *Dv1 = kick->Dv1[l] * u + kick->Dv1[l + 1] * v;
        *Dv2 = kick->Dv2[l] * u + kick->Dv2[l + 1] * v;
    }
}

inline void
fastpm_kick_one(FastPMKickFactor * kick, FastPMStore * p, ptrdiff_t i, float vo[3], double af)
{
    double dda_i, Dv1_i, Dv2_i;
    double dda_f, Dv1_f, Dv2_f;
    double dda, Dv1, Dv2;

    fastpm_kick_lookup(kick, af, &dda_f, &Dv1_f, &Dv2_f);
    fastpm_kick_lookup(kick, p->meta.a_v, &dda_i, &Dv1_i, &Dv2_i);
    dda = dda_f - dda_i;
    Dv1 = Dv1_f - Dv1_i;
    Dv2 = Dv2_f - Dv2_i;

    int d;
    for(d = 0; d < 3; d++) {
        float ax = p->acc[i][d];       // unlike a_x, which means a at which x is calcd
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
    po->meta.a_v = af;
}

static double G_p(FastPMGrowthInfo * growth_info)
{
    /* integral of G_p */
    return growth_info->D1;
}
static double g_p(FastPMGrowthInfo * growth_info)
{
    return DGrowthFactorDa(growth_info);
}

static double G_f(FastPMGrowthInfo * growth_info)
{
    /* integral of g_f */
    double a = growth_info->a;
    return a * a * a * HubbleEa(a, growth_info->c) * g_p(growth_info);
}

static double g_f(FastPMGrowthInfo * growth_info)
{
    double a = growth_info->a;
    FastPMCosmology * c = growth_info->c;

    double E = HubbleEa(a, c);
    double dEda = DHubbleEaDa(a, c);

    double dDda = g_p(growth_info);
    double d2Dda2 = D2GrowthFactorDa2(growth_info);

    double g_f = 3 * a * a * E * dDda
                   + a * a * a * dEda * dDda
                   + a * a * a * E * d2Dda2;
    return g_f;
}

void fastpm_kick_init(FastPMKickFactor * kick, FastPMSolver * fastpm, double ai, double ac, double af)
{
    FastPMCosmology * c = fastpm->cosmology;
    kick->forcemode = fastpm->config->FORCE_TYPE;

    FastPMGrowthInfo gi_i;
    FastPMGrowthInfo gi_c;
    FastPMGrowthInfo gi_e;

    fastpm_growth_info_init(&gi_i, ai, c);
    fastpm_growth_info_init(&gi_c, ac, c);

    double E_i = HubbleEa(ai, c);
    double E_c = HubbleEa(ac, c);

    double D1_i = gi_i.D1;
    double D2_i = gi_i.D2;
    double f1_i = gi_i.f1;
    double f2_i = gi_i.f2;

    double D1_c = gi_c.D1;
    double D2_c = gi_c.D2;

    double Omega_m0 = Omega_source(1, c);
    double Omega_mc = Omega_source(ac, c);

    // kick->q1,2 are used for the COLA force implementation.
    // growth_mode = ODE and LCDM should match for an LCDM background,
    // but neither is guaranteed accurate for a background with radiaiton.
    // We advise using LCDM mode for forcemode = FASTPM_FORCE_COLA, as in the
    // original implementation of FastPM.
    kick->q1 = D1_c;
    switch (c->growth_mode){
        case FASTPM_GROWTH_MODE_LCDM:
            kick->q2 = D1_c*D1_c * (1.0 + 7.0/3.0 * pow(Omega_mc, 1.0/143.0));
        break;
        case FASTPM_GROWTH_MODE_ODE:
            kick->q2 = D1_c*D1_c * (1 - D1_c*D1_c/D2_c);
        break;
        default:
            fastpm_raise(-1, "Please enter a valid growth mode.\n");
    }

    kick->nsamples = 32;
    int i;

    double Dv1i = D1_i * ai * ai * E_i * f1_i;
    double Dv2i = D2_i * ai * ai * E_i * f2_i;
    for(i = 0; i < kick->nsamples; i ++) {
        double ae = ai * (1.0 * (kick->nsamples - 1 - i) / (kick->nsamples - 1))
                  + af * (1.0 * i / (kick->nsamples - 1));

        fastpm_growth_info_init(&gi_e, ae, c);
        double D1_e = gi_e.D1;
        double f1_e = gi_e.f1;
        double D2_e = gi_e.D2;
        double f2_e = gi_e.f2;
        double E_e = HubbleEa(ae, c);

        if(kick->forcemode == FASTPM_FORCE_FASTPM) {
            kick->dda[i] = -1.5 * Omega_mc * ac
               * E_c
               * (G_f(&gi_e) - G_f(&gi_i)) / g_f(&gi_c);
        } else {
            kick->dda[i] = -1.5 * Omega_m0
                * Sphi(ai, ae, ac, fastpm->config->nLPT, c, kick->forcemode == FASTPM_FORCE_COLA);
        }
        kick->Dv1[i] = D1_e * ae * ae * E_e * f1_e - Dv1i;
        kick->Dv2[i] = D2_e * ae * ae * E_e * f2_e - Dv2i;
    }

    kick->ai = ai;
    kick->ac = ac;
    kick->af = af;

    /* Output growth and FastPM factor at af for reference.
    This is a weird place to put this, but it's convenient because G and g are static */
    fastpm_info("Growth/FastPM factors at a = %6.4f: D1=%g, D2=%g, f1=%g, f2=%g, G_p=%g, G_f=%g, g_p=%g, g_f=%g\n",
               ai,
               gi_i.D1,
               gi_i.D2,
               gi_i.f1,
               gi_i.f2,
               G_p(&gi_i),
               G_f(&gi_i),
               g_p(&gi_i),
               g_f(&gi_i));
}

void
fastpm_drift_init(FastPMDriftFactor * drift, FastPMSolver * fastpm,
                double ai, double ac, double af)
{
    FastPMCosmology * c = fastpm->cosmology;
    drift->forcemode = fastpm->config->FORCE_TYPE;

    FastPMGrowthInfo gi_i;
    FastPMGrowthInfo gi_c;
    FastPMGrowthInfo gi_e;

    fastpm_growth_info_init(&gi_i, ai, c);
    fastpm_growth_info_init(&gi_c, ac, c);

    double E_c = HubbleEa(ac, c);

    double D1_i = gi_i.D1;
    double D2_i = gi_i.D2;

    double D1_c = gi_c.D1;
    double D2_c = gi_c.D2;
    double f1_c = gi_c.f1;
    double f2_c = gi_c.f2;

    drift->nsamples = 32;
    int i;

    for(i = 0; i < drift->nsamples; i ++ ) {
        double ae = ai * (1.0 * (drift->nsamples - 1 - i) / (drift->nsamples - 1))
                  + af * (1.0 * i / (drift->nsamples - 1));

        fastpm_growth_info_init(&gi_e, ae, c);     // overwrite each iteration
        double D1_e = gi_e.D1;
        double D2_e = gi_e.D2;

        if (drift->forcemode == FASTPM_FORCE_FASTPM) {
            drift->dyyy[i] = 1 / (ac * ac * ac * E_c)
                        * (G_p(&gi_e) - G_p(&gi_i)) / g_p(&gi_c);
        } else {
            drift->dyyy[i] = Sq(ai, ae, ac, fastpm->config->nLPT, c, drift->forcemode == FASTPM_FORCE_COLA);
        }
        drift->da1[i] = D1_e - D1_i;    // change in D_1lpt
        drift->da2[i] = D2_e - D2_i;  // change in D_2lpt
    }
    drift->af = af;
    drift->ai = ai;
    drift->ac = ac;
    drift->Dv1 = D1_c * ac * ac * E_c * f1_c;
    drift->Dv2 = D2_c * ac * ac * E_c * f2_c;
}

void
fastpm_drift_store(FastPMDriftFactor * drift,
               FastPMStore * pi, FastPMStore * po,
               double af)
{
    int np = pi->np;

    int i;
    // Drift
#pragma omp parallel for
    for(i=0; i<np; i++) {
        double xo[3] = {0};
        fastpm_drift_one(drift, pi, i, xo, af);
        int d;
        for(d = 0; d < 3; d ++) {
            po->x[i][d] = xo[d];
        }
    }
    po->meta.a_x = af;
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
