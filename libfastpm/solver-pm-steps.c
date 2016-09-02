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
#include <fastpm/cosmology.h>

#include "pmpfft.h"
#include "vpm.h"

static double 
Sq(double ai, double af, double aRef, FastPMSolver * );

static double 
Sphi(double ai, double af, double aRef, FastPMSolver * );

static Cosmology CP(FastPMSolver * fastpm) {
    Cosmology c = {
        .OmegaM = fastpm->omega_m,
        .OmegaLambda = 1 - fastpm->omega_m,
    };
    return c;
}

double
fastpm_growth_factor(FastPMSolver * fastpm, double a)
{
    return GrowthFactor(a, CP(fastpm));
}

inline void
fastpm_drift_one(FastPMDrift * drift, ptrdiff_t i, double xo[3], double af)
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
        dyyy = drift->dyyy[l] * u + drift->dyyy[l + 1] * v;
        da1  = drift->da1[l] * u + drift->da1[l + 1] * v;
        da2  = drift->da2[l] * u + drift->da2[l + 1] * v;
    }
    int d;
    for(d = 0; d < 3; d ++) {
        if(drift->fastpm->FORCE_TYPE == FASTPM_FORCE_2LPT) {
            xo[d] = drift->p->x[i][d] + drift->p->dx1[i][d] * da1 + drift->p->dx2[i][d] * da2;
        } else if(drift->fastpm->FORCE_TYPE == FASTPM_FORCE_ZA) {
            xo[d] = drift->p->x[i][d] + drift->p->dx1[i][d] * da1;
        } else {
            xo[d] = drift->p->x[i][d] + drift->p->v[i][d] * dyyy;
            if(drift->fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
                xo[d] += drift->p->dx1[i][d] * da1 + drift->p->dx2[i][d] * da2;
            }
        }
    }
}

inline void
fastpm_kick_one(FastPMKick * kick, ptrdiff_t i, float vo[3])
{
    int d;
    for(d = 0; d < 3; d++) {
        float ax = kick->p->acc[i][d];
        if(kick->fastpm->FORCE_TYPE == FASTPM_FORCE_COLA) {
            ax += (kick->p->dx1[i][d]*kick->q1 + kick->p->dx2[i][d]*kick->q2);
        }
        vo[d] = kick->p->v[i][d] + ax * kick->dda;
    }
}

// Leap frog time integration

void 
fastpm_kick_store(FastPMSolver * fastpm, 
              FastPMStore * pi, FastPMStore * po, double af)
{
    FastPMKick * kick = alloca(sizeof(FastPMKick));

    fastpm_kick_init(kick, fastpm, pi, af);

    int np = pi->np;

    // Kick using acceleration at a= ac
    // Assume forces at a=ac is in particles->force

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        float vo[3];
        fastpm_kick_one(kick, i, vo);
        for(d = 0; d < 3; d++) {
            po->v[i][d] = vo[d];
        }
    }

    //velocity is now at a= avel1
    po->a_v = af;
}

static double G_p(double a, Cosmology c)
{
    /* integral of G_p */
    return GrowthFactor(a, c);
}
static double g_p(double a, Cosmology c)
{
    return DGrowthFactorDa(a, c);
}

static double G_f(double a, Cosmology c)
{
    /* integral of g_f */
    return a * a * a * HubbleEa(a, c) * g_p(a, c);
}

static double g_f(double a, Cosmology c)
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
void fastpm_kick_init(FastPMKick * kick, FastPMSolver * fastpm, FastPMStore * pi, double af)
{
    double ai = pi->a_v;
    double ac = pi->a_x;
    double OmegaM = fastpm->omega_m;
    Cosmology c = CP(fastpm);
    double Om143 = pow(OmegaA(ac, c), 1.0/143.0);
    double growth1 = GrowthFactor(ac, c);

    kick->q2 = growth1*growth1*(1.0 + 7.0/3.0*Om143);
    kick->q1 = growth1;
    if(fastpm->FORCE_TYPE == FASTPM_FORCE_FASTPM) {
        kick->dda = -1.5 * OmegaM
           * 1 / (ac * ac * HubbleEa(ac, c))
           * (G_f(af, c) - G_f(ai, c)) / g_f(ac, c);
    } else {
        kick->dda = -1.5 * OmegaM * Sphi(ai, af, ac, fastpm);
    }
    kick->fastpm = fastpm;
    kick->p = pi;
    kick->af = af;
}

void
fastpm_drift_init(FastPMDrift * drift, FastPMSolver * fastpm, FastPMStore * pi,
               double af)
{
    double ai = pi->a_x;
    double ac = pi->a_v;

    Cosmology c = CP(fastpm);

    drift->nsamples = 32;
    int i;

    for(i = 0; i < drift->nsamples; i ++ ) {
        double ae = ai * (1.0 * (drift->nsamples - 1 - i) / (drift->nsamples - 1))
                  + af * (1.0 * i / (drift->nsamples - 1));

        if(fastpm->FORCE_TYPE == FASTPM_FORCE_FASTPM) {
            drift->dyyy[i] = 1 / (ac * ac * ac * HubbleEa(ac, c))
                        * (G_p(ae, c) - G_p(ai, c)) / g_p(ac, c);
        } else {
            drift->dyyy[i] = Sq(ai, ae, ac, fastpm);
        }
        drift->da1[i] = GrowthFactor(ae, c) - GrowthFactor(ai, c);    // change in D_1lpt
        drift->da2[i] = GrowthFactor2(ae, c) - GrowthFactor2(ai, c);  // change in D_2lpt
    }
    drift->fastpm = fastpm;
    drift->p = pi;
    drift->af = af;
    drift->ai = ai;
}

void
fastpm_drift_store(FastPMSolver * fastpm,
               FastPMStore * pi, FastPMStore * po,
               double af)
{

    FastPMDrift * drift = alloca(sizeof(FastPMDrift));

    fastpm_drift_init(drift, fastpm, pi, af);

    int np = pi->np;

    int i;
    // Drift
#pragma omp parallel for
    for(i=0; i<np; i++) {
        double xo[3];
        fastpm_drift_one(drift, i, xo, drift->af);
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

double gpQ(double a, double nLPT) { 
    return pow(a, nLPT);
}

static double stddriftfunc (double a, FastPMSolver * fastpm) {
    return 1 / (pow(a, 3) * HubbleEa(a, CP(fastpm)));
}

static double nonstddriftfunc (double a, FastPMSolver * fastpm) {
    return gpQ(a, fastpm->nLPT)/(pow(a, 3) * HubbleEa(a, CP(fastpm)));
}

static double stdkickfunc (double a, FastPMSolver * fastpm) {
    return 1/ (pow(a, 2) * HubbleEa(a, CP(fastpm)));
}

static double integrand(double a, void * params) {
    void ** p = (void**) params;
    double (*func)(double a, FastPMSolver * s) = p[0];
    FastPMSolver * s = p[1];
    return func(a, s);
}

double integrate(double ai, double af,
        FastPMSolver * fastpm,
        double (*func)(double , FastPMSolver * )) {

    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (5000);

    gsl_function F;
    double error;
    double result;

    F.params = (void*[]){func, fastpm};
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
Sq(double ai, double af, double aRef, FastPMSolver * fastpm)
{
    double resultstd, result;

    resultstd = integrate(ai, af, fastpm, stddriftfunc);

    result = integrate(ai, af, fastpm, nonstddriftfunc);
    result /= gpQ(aRef, fastpm->nLPT);

    /*
    fastpm_info("ref time = %6.4f, std drift =%g, non std drift = %g \n",
        aRef, resultstd, result); */

    if (fastpm->USE_NONSTDDA)
        return result;
    else
        return resultstd;
}

double DERgpQ(double a, double nLPT) { 
    /* This must return d(gpQ)/da */
    return nLPT*pow(a, nLPT-1);
}


static double 
Sphi(double ai, double af, double aRef, FastPMSolver * fastpm) 
{
    double result;
    double resultstd;

    result = (gpQ(af, fastpm->nLPT) - gpQ(ai, fastpm->nLPT)) * aRef 
        / (pow(aRef, 3) * HubbleEa(aRef, CP(fastpm)) * DERgpQ(aRef, fastpm->nLPT));

    resultstd = integrate(ai, af, fastpm, stdkickfunc);

    /*
    fastpm_info("ref time = %6.4f, std kick = %g, non std kick = %g\n",
            aRef, resultstd, result); */

    if (fastpm->USE_NONSTDDA) {
        return result;
    } else {
        return resultstd;
    }
}


// Interpolate position and velocity for snapshot at a=aout
void 
fastpm_set_snapshot(FastPMSolver * fastpm,
                FastPMStore * p, FastPMStore * po,
                double aout)
{
    int np= p->np;
    double a_x = p->a_x;
    double a_v = p->a_v;

    fastpm_info("Setting up snapshot at a= %6.4f (z=%6.4f) <- %6.4f %6.4f.\n", aout, 1.0f/aout-1, a_x, a_v);

    double H0 = 100.0f; // H0= 100 km/s/(h^-1 Mpc)

    fastpm_info("Growth factor of snapshot %f (a=%.3f)\n", fastpm_growth_factor(fastpm, aout), aout);

    Cosmology c = CP(fastpm);

    double Dv1 = GrowthFactor(aout, c) * aout * aout * HubbleEa(aout, c) * DLogGrowthFactor(aout, c);
    double Dv2 = GrowthFactor2(aout, c) * aout * aout * HubbleEa(aout, c) * DLogGrowthFactor2(aout, c);

    fastpm_info("RSD factor %e\n", 1 /(aout * HubbleEa(aout, c) *H0));

    fastpm_kick_store(fastpm, p, po, aout);

    fastpm_drift_store(fastpm, p, po, aout);

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* For cola,
             * add the lpt velocity to the residual velocity v*/
            if(fastpm->FORCE_TYPE == FASTPM_FORCE_COLA)
                po->v[i][d] += p->dx1[i][d]*Dv1
                             + p->dx2[i][d]*Dv2;
            /* convert the unit from a**2 H_0 dx/dt in Mpc/h to a dx/dt km/s */
            po->v[i][d] *= H0 / aout;
        }
        po->id[i] = p->id[i];
    }

    po->np = np;
    po->a_x = po->a_v = aout;
}

