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

#include "pmpfft.h"
#include "parameters.h"
#include "msg.h"
#include "pmsteps.h"
#include "cosmology.h"

static double nLPT= -2.5f;
static int martinKick = 0;

static double stepping_boost = 1.0;

static double 
Sq(double ai, double af, double aRef, PMStepper * );

static double 
Sphi(double ai, double af, double aRef, PMStepper * );

static Cosmology stepper_get_cosmology(PMStepper * stepper) {
    Cosmology c = {
        .OmegaM = stepper->omega_m,
        .OmegaLambda = 1 - stepper->omega_m,
    };
    return c;
}
double stepper_get_growth_factor(PMStepper * stepper, double a) {
    Cosmology c = {
        .OmegaM = stepper->omega_m,
        .OmegaLambda = 1 - stepper->omega_m,
    };
    double growth1 = GrowthFactor(a, c);
    return growth1;
}

void stepping_set_boost(double boost) {
    stepping_boost = boost;
}

void 
stepping_init(PMStepper * stepper, double omega_m, int force_mode, int stdDA) 
{
    stepper->omega_m = omega_m;
    stepper->mode = force_mode;
    if(force_mode == FORCE_MODE_COLA) {
        stepper->stdda = stdDA;
    } else {
        stepper->stdda = 1;
    }
}

// Leap frog time integration
// ** Total momentum adjustment dropped

void 
stepping_kick(PMStepper * stepper, 
              PMStore * pi, PMStore * po,
              double ai, double af, double ac)
                /* a_v     avel1     a_x*/
{
    Cosmology c = stepper_get_cosmology(stepper);
    double Om143 = pow(OmegaA(ac, c), 1.0/143.0);
    double dda = Sphi(ai, af, ac, stepper) * stepping_boost;
    double growth1 = GrowthFactor(ac, c);

    msg_printf(normal, "Kick %6.4f -> %6.4f\n", ai, af);
    msg_printf(normal, "growth factor %g dda=%g\n", growth1, dda);

    double q2 = 1.5*c.OmegaM*growth1*growth1*(1.0 + 7.0/3.0*Om143);
    double q1 = 1.5*c.OmegaM*growth1;

    double omegam;
    if (martinKick)
        omegam = OmegaA(ac, c);
    else
        omegam = c.OmegaM;

    int np = pi->np;

    // Kick using acceleration at a= ac
    // Assume forces at a=ac is in particles->force

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d++) {
            float ax= -1.5 * omegam * pi->acc[i][d];
            switch(stepper->mode) {
                case FORCE_MODE_COLA:
                    ax -= (pi->dx1[i][d]*q1 + pi->dx2[i][d]*q2);
                break;
            }
            po->v[i][d] = pi->v[i][d] + ax * dda;
        }
    }

    //velocity is now at a= avel1
}

void 
stepping_drift(PMStepper * stepper,
               PMStore * pi, PMStore * po,
               double ai, double af, double ac)
               /*a_x, apos1, a_v */
{
    int np = pi->np;

    Cosmology c = stepper_get_cosmology(stepper);

    double dyyy = Sq(ai, af, ac, stepper) * stepping_boost;

    double da1 = GrowthFactor(af, c) - GrowthFactor(ai, c);    // change in D_1lpt
    double da2 = GrowthFactor2(af, c) - GrowthFactor2(ai, c);  // change in D_2lpt

    msg_printf(normal, "Drift %6.4f -> %6.4f\n", ai, af);
    msg_printf(normal, "dyyy = %g \n", dyyy);


    int i;
    // Drift
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            po->x[i][d] = pi->x[i][d] + pi->v[i][d]*dyyy;
            switch(stepper->mode) {
                case FORCE_MODE_COLA:
                    po->x[i][d] += pi->dx1[i][d]*da1 + pi->dx2[i][d]*da2;
                break;
            }
        }
    }

}

//
// Functions for our modified time-stepping (used when StdDA=0):
//

double gpQ(double a) { 
    return pow(a, nLPT);
}

static double stddriftfunc (double a, PMStepper * stepper) {
    return 1.0/Qfactor(a, stepper_get_cosmology(stepper));
}

static double nonstddriftfunc (double a, PMStepper * stepper) {
    return gpQ(a)/Qfactor(a, stepper_get_cosmology(stepper)); 
}

static double stdkickfunc (double a, PMStepper * stepper) {
    return a/Qfactor(a, stepper_get_cosmology(stepper));
}

static double integrand(double a, void * params) {
    void ** p = (void**) params;
    double (*func)(double a, PMStepper * s) = p[0];
    PMStepper * s = p[1];
    return func(a, s);
}

double integrate(double ai, double af,
        PMStepper * stepper,
        double (*func)(double , PMStepper * )) {

    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (5000);
    
    gsl_function F;
    double error;
    double result;

    F.params = (void*[]){func, stepper};
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
Sq(double ai, double af, double aRef, PMStepper * stepper)
{
    double resultstd, result;

    resultstd = integrate(ai, af, stepper, stddriftfunc);

    result = integrate(ai, af, stepper, nonstddriftfunc);
    result /= gpQ(aRef);

    msg_printf(verbose, "ref time = %6.4f, std drift =%g, non std drift = %g \n",
        aRef, resultstd, result);

    if (stepper->stdda == 0)
        return result;
    else
        return resultstd;
}

double DERgpQ(double a) { // This must return d(gpQ)/da
    return nLPT*pow(a, nLPT-1);
}


static double 
Sphi(double ai, double af, double aRef, PMStepper * stepper) 
{
    double result;
    double resultstd;

    Cosmology c = stepper_get_cosmology(stepper);

    /* Qfactor is a**2 da / dt */
    result = (gpQ(af) - gpQ(ai)) * aRef 
        / (Qfactor(aRef, c) * DERgpQ(aRef));

    resultstd = integrate(ai, af, stepper, stdkickfunc);

    msg_printf(verbose, "ref time = %6.4f, std kick = %g, non std kick = %g\n",
            aRef, resultstd, result);

    if (stepper->stdda == 0) {
        return result;
    } else {
        return resultstd;
    }
}


// Interpolate position and velocity for snapshot at a=aout
void 
stepping_set_snapshot(PMStepper * stepper,
                PMStore * p, PMStore * po,
                double aout, double a_x, double a_v)
{
    int np= p->np;

    Cosmology c = stepper_get_cosmology(stepper);

    msg_printf(verbose, "Setting up snapshot at a= %6.4f (z=%6.4f) <- %6.4f %6.4f.\n", aout, 1.0f/aout-1, a_x, a_v);

    float vfac= 100.0f/aout;   // km/s; H0= 100 km/s/(h^-1 Mpc)

    msg_printf(normal, "Growth factor of snapshot %f (a=%.3f)\n", GrowthFactor(aout, c), aout);

    double Dv=DprimeQ(aout, 1.0, c); // dD_{za}/dy
    double Dv2=GrowthFactor2v(aout, c);   // dD_{2lpt}/dy

    msg_printf(debug, "velocity factor %e %e\n", vfac*Dv, vfac*Dv2);
    msg_printf(debug, "RSD factor %e\n", aout/Qfactor(aout, c)/vfac);

    stepping_kick(stepper, p, po, a_v, aout, a_x);

    stepping_drift(stepper, p, po, a_x, aout, a_v);

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* For cola, 
             * add the lpt velocity to the residual velocity v*/
            switch(stepper->mode) {
                case FORCE_MODE_COLA:
                    po->v[i][d] += p->dx1[i][d]*Dv 
                                 + p->dx2[i][d]*Dv2;
                break;
            }
            /* convert the unit to km/s */
            po->v[i][d] *= vfac;
        }    
        po->id[i] = p->id[i];
    }

    po->np = np;
}

