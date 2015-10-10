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
static int stdDA = 0; // velocity growth model
static int martinKick = 0;

double Sphi(double ai, double af, double aRef, Cosmology c);
double Sq(double ai, double af, double aRef, Cosmology c);
static int FORCE_MODE;
static int NSTEPS;
static double *A_X;
static double *A_V;

static double stepping_boost = 1.0;

static double polval(const double * pol, const int degrees, const double x) {
    double rt = 0;
    int i;
    for(i = 0; i < degrees; i ++) {
        rt += pow(x, degrees - i - 1) * pol[i];
    }
    return rt;
}

void stepping_set_boost(double boost) {
    stepping_boost = boost;
}

void stepping_init(Parameters * param) {
    FORCE_MODE = param->force_mode;
    NSTEPS = param->n_time_step;
    stdDA = param->cola_stdda;

    if(param->force_mode == FORCE_MODE_2LPT ||
       param->force_mode == FORCE_MODE_ZA) {
        NSTEPS = 1;
        A_X = (double*) calloc(NSTEPS + 1, sizeof(double));
        A_V = (double*) calloc(NSTEPS + 1, sizeof(double));
        A_V[0] = param->time_step[0];
        A_X[0] = param->time_step[0];
        A_X[NSTEPS] = 1.0;
        A_V[NSTEPS] = 1.0;
    } else {
        /* one extra item in the end; to avoid an if conditioni in main loop */
        A_X = (double*) calloc(NSTEPS + 1, sizeof(double));
        A_V = (double*) calloc(NSTEPS + 1, sizeof(double));
        int i;
        for (i = 0; i <= NSTEPS - 1;i++){
            A_X[i] = param->time_step[i];
        }
        A_X[NSTEPS] = 1.0;

        A_V[0] = A_X[0];

        for (i = 1; i <= NSTEPS; i++){
            A_V[i] = exp((log(A_X[i])+log(A_X[i-1]))/2);
        }
    }
    char * buf = malloc(24 * (NSTEPS + 1));
    char * p;
    int i;
    msg_printf(normal, "%d steps: \n", NSTEPS);
    msg_printf(normal, "Drift: \n");
    for(p = buf, i = 0; i < NSTEPS + 1; i ++) {
        sprintf(p, "%6.4f ", A_X[i]);
        p += strlen(p);
    }
    msg_printf(normal, "%s\n", buf);
    msg_printf(normal, " Kick: \n");
    for(p = buf, i = 0; i < NSTEPS + 1; i ++) {
        sprintf(p, "%6.4f ", A_V[i]);
        p += strlen(p);
    }
    msg_printf(normal, "%s\n", buf);
}

int stepping_get_nsteps() {
    return NSTEPS;
}
void stepping_get_times(int istep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1) {

    *a_v = A_V[istep];
    *a_x = A_X[istep];
    if(istep >= NSTEPS) {
        istep = NSTEPS - 1;
    }
    *a_v1= A_V[istep + 1];
    *a_x1= A_X[istep + 1];
}
// Leap frog time integration
// ** Total momentum adjustment dropped

void 
stepping_kick(PMStore * pi, PMStore * po,
              double ai, double af, double ac,
                /* a_v     avel1     a_x*/
              double OmegaM)
{
    Cosmology c = {
        .OmegaM = OmegaM,
        .OmegaLambda = 1 - OmegaM,
    };
    if(FORCE_MODE == FORCE_MODE_ZA
    || FORCE_MODE == FORCE_MODE_2LPT) {
        /* ZA and 2LPT sims no kicks */
        return;
    }

    double Om143 = pow(OmegaA(ac, c), 1.0/143.0);
    double dda = Sphi(ai, af, ac, c) * stepping_boost;
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
            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                    ax -= (pi->dx1[i][d]*q1 + pi->dx2[i][d]*q2);
                break;
                case FORCE_MODE_COLA1:
                    ax -= (pi->dx1[i][d]*q1);
                break;
            }
            po->v[i][d] = pi->v[i][d] + ax * dda;
        }
    }

    //velocity is now at a= avel1
}

void 
stepping_drift(PMStore * pi, PMStore * po,
               double ai, double af, double ac,
              /*a_x, apos1, a_v */
               double OmegaM)
{
    int np = pi->np;
    Cosmology c = {
        .OmegaM = OmegaM,
        .OmegaLambda = 1 - OmegaM,
    };

    double dyyy = Sq(ai, af, ac, c) * stepping_boost;

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
            if(FORCE_MODE & FORCE_MODE_PM) {
                po->x[i][d] = pi->x[i][d] + pi->v[i][d]*dyyy;
            } else {
                po->x[i][d] = 0;
            }
            switch(FORCE_MODE) {
                case FORCE_MODE_2LPT:
                case FORCE_MODE_COLA:
                    po->x[i][d] += pi->dx1[i][d]*da1 + pi->dx2[i][d]*da2;
                break;
                case FORCE_MODE_ZA:
                case FORCE_MODE_COLA1:
                    po->x[i][d] += pi->dx1[i][d]*da1;
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

double stddriftfunc (double a, void * params) {
    Cosmology * c = params;
    return 1.0/Qfactor(a, *c);
}

double nonstddriftfunc (double a, void * params) {
    Cosmology * c = params;
    return gpQ(a)/Qfactor(a, *c); 
}

double stdkickfunc (double a, void * params) {
    Cosmology * c = params;
    return a/Qfactor(a, *c);
}

double martinkickfunc (double a, void * params) {
    Cosmology * c = params;
    return Qfactor(a, *c) / (a*a);
}

/*     
       When StdDA=0, one needs to set nLPT.
       assumes time dep. for velocity = B a^nLPT
       nLPT is a real number. Sane values lie in the range (-4,3.5). Cannot be 0, but of course can be -> 0 (say 0.001).
       See Section A.3 of TZE.
       */

double Sq(double ai, double af, double aRef, Cosmology c) {
    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (5000);

    double resultstd, result, error;

    gsl_function F;
    F.function = &stddriftfunc;
    F.params = &c;
    gsl_integration_qag (&F, ai, af, 0, 1e-8, 5000,6,
            w, &resultstd, &error); 

    F.function = &nonstddriftfunc;
    F.params = &c;
    gsl_integration_qag (&F, ai, af, 0, 1e-8, 5000,6,
            w, &result, &error); 

    result /= gpQ(aRef);

    gsl_integration_workspace_free (w);

    msg_printf(verbose, "ref time = %6.4f, std drift =%g, non std drift = %g \n",
        aRef, resultstd, result);

    if (stdDA == 0)
        return result;
    else
        return resultstd;
}

double DERgpQ(double a) { // This must return d(gpQ)/da
    return nLPT*pow(a, nLPT-1);
}

double Sphi(double ai, double af, double aRef, Cosmology c) {
    double result;
    double resultstd;

    /* Qfactor is a**2 da / dt */
    result = (gpQ(af) - gpQ(ai)) * aRef 
        / (Qfactor(aRef, c) * DERgpQ(aRef));

    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (5000);

    double error;

    gsl_function F;
    if (martinKick) {
        F.function = &martinkickfunc;
    } else{
        F.function = &stdkickfunc;
    }
    F.params = &c;

    gsl_integration_qag (&F, ai, af, 0, 1e-8, 5000, 6,
            w, &resultstd, &error); 

    gsl_integration_workspace_free (w);

    msg_printf(verbose, "ref time = %6.4f, std kick = %g, non std kick = %g\n",
            aRef, resultstd, result);

    if (stdDA == 0) {
        return result;
    } else {
        return resultstd;
    }
}


// Interpolate position and velocity for snapshot at a=aout
void 
stepping_set_snapshot(PMStore * p, PMStore * po,
                double aout, double a_x, double a_v,
                double OmegaM)
{
    int np= p->np;

    Cosmology c = {
        .OmegaM = OmegaM,
        .OmegaLambda = 1 - OmegaM,
    };

    msg_printf(verbose, "Setting up snapshot at a= %6.4f (z=%6.4f) <- %6.4f %6.4f.\n", aout, 1.0f/aout-1, a_x, a_v);

    float vfac= 100.0f/aout;   // km/s; H0= 100 km/s/(h^-1 Mpc)

    msg_printf(normal, "Growth factor of snapshot %f (a=%.3f)\n", GrowthFactor(aout, c), aout);

    double Dv=DprimeQ(aout, 1.0, c); // dD_{za}/dy
    double Dv2=GrowthFactor2v(aout, c);   // dD_{2lpt}/dy

    msg_printf(debug, "velocity factor %e %e\n", vfac*Dv, vfac*Dv2);
    msg_printf(debug, "RSD factor %e\n", aout/Qfactor(aout, c)/vfac);

    stepping_kick(p, po, a_v, aout, a_x, OmegaM);

    stepping_drift(p, po, a_x, aout, a_v, OmegaM);

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* For cola, 
             * add the lpt velocity to the residual velocity v*/
            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                    po->v[i][d] += p->dx1[i][d]*Dv 
                                 + p->dx2[i][d]*Dv2;
                break;
                case FORCE_MODE_COLA1:
                    po->v[i][d] += p->dx1[i][d]*Dv;
                break;
            }
            /* convert the unit to km/s */
            po->v[i][d] *= vfac;
        }    
        po->id[i] = p->id[i];
    }

    po->np = np;
}

