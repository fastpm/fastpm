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
#include <mpi.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "vpm.h"

static double 
Sq(double ai, double af, double aRef, FastPM * );

static double 
Sphi(double ai, double af, double aRef, FastPM * );

static Cosmology CP(FastPM * fastpm) {
    Cosmology c = {
        .OmegaM = fastpm->omega_m,
        .OmegaLambda = 1 - fastpm->omega_m,
    };
    return c;
}
double 
fastpm_growth_factor(FastPM * fastpm, double a) 
{
    return GrowthFactor(a, CP(fastpm));
}

// Leap frog time integration

void 
fastpm_kick_store(FastPM * fastpm, 
              PMStore * pi, PMStore * po, double af)
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

void fastpm_kick_init(FastPMKick * kick, FastPM * fastpm, PMStore * pi, double af)
{
    double ai = pi->a_v;
    double ac = pi->a_x;
    double OmegaM = fastpm->omega_m;

    double Om143 = pow(OmegaA(ac, CP(fastpm)), 1.0/143.0);
    double growth1 = GrowthFactor(ac, CP(fastpm));

    kick->q2 = growth1*growth1*(1.0 + 7.0/3.0*Om143);
    kick->q1 = growth1;
    kick->dda = -1.5 * OmegaM * Sphi(ai, af, ac, fastpm);
    kick->fastpm = fastpm;
    kick->p = pi;
    kick->af = af;
}

void
fastpm_drift_init(FastPMDrift * drift, FastPM * fastpm, PMStore * pi,
               double af)
{
    double ai = pi->a_x;
    double ac = pi->a_v;

    drift->dyyy = Sq(ai, af, ac, fastpm);

    drift->da1 = GrowthFactor(af, CP(fastpm)) - GrowthFactor(ai, CP(fastpm));    // change in D_1lpt
    drift->da2 = GrowthFactor2(af, CP(fastpm)) - GrowthFactor2(ai, CP(fastpm));  // change in D_2lpt

    drift->fastpm = fastpm;
    drift->p = pi;
    drift->af = af;
}

void
fastpm_drift_store(FastPM * fastpm,
               PMStore * pi, PMStore * po,
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
        fastpm_drift_one(drift, i, xo);
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

static double stddriftfunc (double a, FastPM * fastpm) {
    return 1 / (pow(a, 3) * HubbleEa(a, CP(fastpm)));
}

static double nonstddriftfunc (double a, FastPM * fastpm) {
    return gpQ(a, fastpm->nLPT)/(pow(a, 3) * HubbleEa(a, CP(fastpm)));
}

static double stdkickfunc (double a, FastPM * fastpm) {
    return 1/ (pow(a, 2) * HubbleEa(a, CP(fastpm)));
}

static double integrand(double a, void * params) {
    void ** p = (void**) params;
    double (*func)(double a, FastPM * s) = p[0];
    FastPM * s = p[1];
    return func(a, s);
}

double integrate(double ai, double af,
        FastPM * fastpm,
        double (*func)(double , FastPM * )) {

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
Sq(double ai, double af, double aRef, FastPM * fastpm)
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
Sphi(double ai, double af, double aRef, FastPM * fastpm) 
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
fastpm_set_snapshot(FastPM * fastpm,
                PMStore * p, PMStore * po,
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
            if(fastpm->USE_COLA)
                po->v[i][d] += p->dx1[i][d]*Dv1
                             + p->dx2[i][d]*Dv2;
            /* convert the unit from a**2 dx/dt in Mpc/h to a dx/dt km/s */
            po->v[i][d] *= H0 / aout;
        }
        po->id[i] = p->id[i];
    }

    po->np = np;
    po->a_x = po->a_v = aout;
}

