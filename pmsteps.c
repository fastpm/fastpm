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

static double Omega= -1.0f;
static double OmegaLambda= -1.0f;
static double nLPT= -2.5f;
static int stdDA = 0; // velocity growth model
static int martinKick = 0;

double growthD(double a);
double growthD2(double a);
double Sphi(double ai, double af, double aRef);
double Sq(double ai, double af, double aRef);
double Qfactor(double a);
double GrowthFactor(double astart, double aend);
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
    Omega = param->omega_m;
    OmegaLambda = 1.0 - param->omega_m;

//    double a_final= param->a_final;
//    double a_init = param->a_init;

    if(param->force_mode == FORCE_MODE_2LPT ||
       param->force_mode == FORCE_MODE_ZA) {
        NSTEPS = 1;
        A_X = (double*) calloc(NSTEPS + 2, sizeof(double));
        A_V = (double*) calloc(NSTEPS + 2, sizeof(double));
        A_V[1] = param->time_step[param->n_time_step-1];
        A_X[1] = param->time_step[param->n_time_step-1];
    } else {
        /* one extra item in the end; to avoid an if conditioni in main loop */
        A_X = (double*) calloc(NSTEPS + 2, sizeof(double));
        A_V = (double*) calloc(NSTEPS + 2, sizeof(double));
        int i;
        for (i = 0;i<=param->n_time_step-1;++i){
            A_X[i] = param->time_step[i];
        }

        A_V[0] = A_X[0];

        for (i = 1;i<=param->n_time_step-1;++i){
            A_V[i] = exp((log(A_X[i])+log(A_X[i-1]))/2);
        }
    }
    A_X[NSTEPS+1] = 1.0;
    A_X[NSTEPS] = 1.0;
    A_V[NSTEPS+1] = 1.0;
    A_V[NSTEPS] = 1.0;
    char * buf = malloc(24 * (NSTEPS + 2));
    char * p;
    int i;
    msg_printf(normal, "Drift: \n");
    for(p = buf, i = 0; i < NSTEPS + 2; i ++) {
        sprintf(p, "%6.4f ", A_X[i]);
        p += strlen(p);
    }
    msg_printf(normal, "%s\n", buf);
    msg_printf(normal, " Kick: \n");
    for(p = buf, i = 0; i < NSTEPS + 2; i ++) {
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
    *a_v1= A_V[istep + 1];
    *a_x1= A_X[istep + 1];
}
// Leap frog time integration
// ** Total momentum adjustment dropped

void stepping_kick(PMStore * p,
        double ai, double af, double ac)
    /* a_v     avel1     a_x*/
{
    if(FORCE_MODE == FORCE_MODE_ZA
    || FORCE_MODE == FORCE_MODE_2LPT) {
        /* ZA and 2LPT sims no kicks */
        return;
    }
    msg_printf(normal, "Kick %6.4f -> %6.4f\n", ai, af);

    double Om143= pow(Omega/(Omega + (1 - Omega)*ac*ac*ac), 1.0/143.0);
    double dda= Sphi(ai, af, ac) * stepping_boost;
    double growth1=growthD(ac);

    msg_printf(normal, "growth factor %g dda=%g\n", growth1, dda);

    double q2=1.5*Omega*growth1*growth1*(1.0 + 7.0/3.0*Om143);
    double q1=1.5*Omega*growth1;

    double Om_;
    if (martinKick)
        Om_ = Omega/ (Omega + (1 - Omega) *ac * ac *ac);
    else
        Om_ = Omega;

    int np = p->np;

    // Kick using acceleration at a= ac
    // Assume forces at a=ac is in particles->force

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d++) {
            float ax= -1.5*Om_* p->acc[i][d];
            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                    ax -= (p->dx1[i][d]*q1 + p->dx2[i][d]*q2);
                break;
                case FORCE_MODE_COLA1:
                    ax -= (p->dx1[i][d]*q1);
                break;
            }
            p->v[i][d] += ax * dda;
        }
    }

    //velocity is now at a= avel1
}

void stepping_drift(PMStore * p,
        double ai, double af, double ac)
    /*a_x, apos1, a_v */
{
    int np= p->np;


    double dyyy=Sq(ai, af, ac);

    double da1= growthD(af) - growthD(ai);    // change in D_1lpt
    double da2= growthD2(af) - growthD2(ai);  // change in D_2lpt

    msg_printf(normal, "Drift %6.4f -> %6.4f\n", ai, af);
    msg_printf(normal, "dyyy = %g \n", dyyy);

    dyyy *= stepping_boost;

    int i;
    // Drift
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            if(FORCE_MODE & FORCE_MODE_PM) {
                p->x[i][d] += p->v[i][d]*dyyy;
            }
            switch(FORCE_MODE) {
                case FORCE_MODE_2LPT:
                case FORCE_MODE_COLA:
                    p->x[i][d] += p->dx1[i][d]*da1 + p->dx2[i][d]*da2;
                break;
                case FORCE_MODE_ZA:
                case FORCE_MODE_COLA1:
                    p->x[i][d] += p->dx1[i][d]*da1;
                break;
            }
        }
    }

}

double growthDtemp(double a){
    // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration?
    double x=-Omega/(Omega - 1.0)/(a*a*a);


    double hyperP=0,hyperM=0;

    if (fabs(x-1.0) < 1.e-3) {
        hyperP= 0.859596768064608 - 0.1016599912520404*(-1.0 + x) + 0.025791094277821357*pow(-1.0 + x,2) - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) - 0.0001788889964687831*pow(-1.0 + x,7);
        hyperM= 1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) +    0.000025063737576245116*pow(-1.0 + x,7);
    }
    else {
        if (x < 1.0) {
            hyperP=gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
            hyperM=gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
        }
        x=1.0/x;
        if ((x < 1.0) && (x>1.0/30)) {

            hyperP=gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
            hyperP*=4*sqrt(x);
            hyperP+=-3.4494794123063873799*pow(x,2.0/3.0);

            hyperM=gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
            hyperM*=4.0/7.0/sqrt(x);
            hyperM+=pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
        }
        if (x<=1.0/30.0){
            hyperP=3.9999999999999996*sqrt(x) - 3.4494794123063865*pow(x,0.6666666666666666) + 0.3999999999999999*pow(x,1.5) -    0.13636363636363635*pow(x,2.5) + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) +    0.033943965517241374*pow(x,5.5) - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) -    0.01671324384973404*pow(x,8.5) + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) - 0.009080577904069926*pow(x,12.5);
            hyperM=0.5714285714285715/sqrt(x) + 2.000000000000001*sqrt(x) - 1.4783483195598794*pow(x,0.66666666666666666) +    0.10000000000000002*pow(x,1.5) - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) -    0.004755434782608697*pow(x,4.5) + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) +    0.0012772722942073172*pow(x,7.5) - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) -    0.0005429801294359111*pow(x,10.5) + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
        }
    }


    if (a > 0.2) 
        return sqrt(1.0 + (-1.0 + pow(a,-3))*Omega)*(3.4494794123063873799*pow(-1.0 + 1.0/Omega,0.666666666666666666666666666) + (hyperP*(4*pow(a,3)*(-1.0 + Omega) - Omega) - 7.0*pow(a,3)*hyperM*(-1.0 + Omega))/(pow(a,5)*(-1.0+ Omega) - pow(a,2)*Omega));

    return (a*pow(1 - Omega,1.5)*(1291467969*pow(a,12)*pow(-1 + Omega,4) + 1956769650*pow(a,9)*pow(-1 + Omega,3)*Omega + 8000000000*pow(a,3)*(-1 + Omega)*pow(Omega,3) + 37490640625*pow(Omega,4)))/(1.5625e10*pow(Omega,5));    
}

double growthD(double a) { // growth factor for LCDM
    return growthDtemp(a)/growthDtemp(1.0);
}


double Qfactor(double a) { // Q\equiv a^3 H(a)/H0.
    return sqrt(Omega/(a*a*a)+1.0-Omega)*a*a*a;
}




double growthD2temp(double a){
    double d= growthD(a);
    double omega=Omega/(Omega + (1.0 - Omega)*a*a*a);
    return d*d*pow(omega, -1.0/143.);
}

double growthD2(double a) {// Second order growth factor
    return growthD2temp(a)/growthD2temp(1.0); // **???
}


double growthD2v(double a){ // explanation is in main()
    double d2= growthD2(a);
    double omega=Omega/(Omega + (1.0 - Omega)*a*a*a);
    return Qfactor(a)*(d2/a)*2.0*pow(omega, 6.0/11.);
}

double decayD(double a){ // D_{-}, the decaying mode
    return sqrt(Omega/(a*a*a)+1.0-Omega);
}

double DprimeQ(double a,double nGrowth)
{ // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
    double nDecay=0.0;// not interested in decay modes in this code.
    double Nn=6.0*pow(1.0 - Omega,1.5)/growthDtemp(1.0);
    return (pow(decayD(a),-1.0 + nDecay)*pow(growthD(a),-1.0 + nGrowth)*(nGrowth*Nn- (3.0*(nDecay + nGrowth)*Omega*growthD(a))/(2.*a)));  
}



//
// Functions for our modified time-stepping (used when StdDA=0):
//

double gpQ(double a) { 
    return pow(a, nLPT);
}

double stddriftfunc (double a, void * params) {
    return 1.0/Qfactor(a);
}

double nonstddriftfunc (double a, void * params) {
    return gpQ(a)/Qfactor(a); 
}

double stdkickfunc (double a, void * params) {
    return a/Qfactor(a);
}

double martinkickfunc (double a, void * params) {
    return Qfactor(a) / (a*a);
}

/*     
       When StdDA=0, one needs to set nLPT.
       assumes time dep. for velocity = B a^nLPT
       nLPT is a real number. Sane values lie in the range (-4,3.5). Cannot be 0, but of course can be -> 0 (say 0.001).
       See Section A.3 of TZE.
       */

double Sq(double ai, double af, double aRef) {
    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (5000);

    double resultstd, result2, error;
    double alpha=0;

    gsl_function F;
    F.function = &stddriftfunc;
    F.params = &alpha;

    gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
            w, &resultstd, &error); 

    F.function = &nonstddriftfunc;
    F.params = &alpha;

    gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
            w, &result2, &error); 
    result2 /= gpQ(aRef);

    gsl_integration_workspace_free (w);

    msg_printf(verbose, "time = %6.4f, std drift =%g, non std drift = %g \n",
        aRef, resultstd, result2);

    if (stdDA == 0)
        return result2;
    else
        return resultstd;
}

double DERgpQ(double a) { // This must return d(gpQ)/da
    return nLPT*pow(a, nLPT-1);
}

double Sphi(double ai, double af, double aRef) {
    double result;
    double resultstd;

    /* Qfactor is a**2 da / dt */
    result = (gpQ(af) - gpQ(ai)) * aRef 
        / (Qfactor(aRef) * DERgpQ(aRef));

    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (5000);

    double error;
    double alpha=0;

    gsl_function F;
    if (martinKick) {
        F.function = &martinkickfunc;
    } else{
        F.function = &stdkickfunc;
    }
    F.params = &alpha;

    gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
            w, &resultstd, &error); 

    gsl_integration_workspace_free (w);

    msg_printf(verbose, "time = %6.4f, std kick = %g, non std kick = %g\n",
            aRef, resultstd, result);

    if (stdDA == 0) {
        return result;
    } else {
        return resultstd;
    }
}


// Interpolate position and velocity for snapshot at a=aout
void stepping_set_initial(double aout, PMStore * p, double shift[3])
{
    int np= p->np;

    msg_printf(verbose, "Setting up inital snapshot at a= %6.4f (z=%6.4f).\n", aout, 1.0f/aout-1);

    const float Dplus = 1.0/GrowthFactor(aout, 1.0);

    const double omega=Omega/(Omega + (1.0 - Omega)*aout*aout*aout);
    const double D2 = Dplus*Dplus*pow(omega/Omega, -1.0/143.0);
    const double D20 = pow(Omega, -1.0/143.0);
    

    float Dv=DprimeQ(aout, 1.0); // dD_{za}/dy
    float Dv2=growthD2v(aout);   // dD_{2lpt}/dy

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* Use the more accurate 2LPT dx2 term */
            p->dx2[i][d] *= 3.0 / 7.0 * D20;

            p->x[i][d] += Dplus * p->dx1[i][d] + D2 * p->dx2[i][d];
            p->x[i][d] += shift[d];

            if(FORCE_MODE != FORCE_MODE_PM) {
                // for COLA and COLA1, v cancels out such that the initial
                // is zero
                p->v[i][d] = 0;
            } else {
                p->v[i][d] = (p->dx1[i][d]*Dv + p->dx2[i][d]*Dv2);
            }
        }
    }
}

// Interpolate position and velocity for snapshot at a=aout
void stepping_set_snapshot(double aout, double a_x, double a_v, 
        PMStore * p, PMStore * po)
{
    int np= p->np;

    msg_printf(verbose, "Setting up snapshot at a= %6.4f (z=%6.4f) <- %6.4f %6.4f.\n", aout, 1.0f/aout-1, a_x, a_v);

    //float vfac=A/Qfactor(A); // RSD /h Mpc unit
    float vfac= 100.0f/aout;   // km/s; H0= 100 km/s/(h^-1 Mpc)

    float AI=  a_v;
    float A=   a_x;
    float AF=  aout;

    float Om143= pow(Omega/(Omega + (1 - Omega)*A*A*A), 1.0/143.0);
    float dda= Sphi(AI, AF, A) * stepping_boost;
    float growth1=growthD(A);

    //msg_printf(normal, "set snapshot %f from %f %f\n", aout, AI, A);
    //msg_printf(normal, "Growth factor of snapshot %f (a=%.3f)\n", growth1, A);
    msg_printf(normal, "Growth factor of snapshot %f (a=%.3f)\n", growthD(AF), AF);

    float q1=1.5*Omega*growth1;
    float q2=1.5*Omega*growth1*growth1*(1.0 + 7.0/3.0*Om143);

    float Dv=DprimeQ(aout, 1.0); // dD_{za}/dy
    float Dv2=growthD2v(aout);   // dD_{2lpt}/dy


    float AC= a_v;
    float dyyy=Sq(A, AF, AC);

    dyyy *= stepping_boost;

    /*
       if(AF < A) {
       float dyyy_backward= Sq(aminus, aout, AC) - Sq(aminus, A, AC);
       msg_printf(debug, "dyyy drift backward %f ->%f = %e %e\n", 
       A, AF, dyyy, dyyy_backward);
       }
       */

    msg_printf(debug, "velocity factor %e %e\n", vfac*Dv, vfac*Dv2);
    msg_printf(debug, "RSD factor %e\n", aout/Qfactor(aout)/vfac);


    float da1= growthD(AF) - growthD(A);    // change in D_{1lpt}
    float da2= growthD2(AF) - growthD2(A);  // change in D_{2lpt}

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            // Kick + adding back 2LPT velocity + convert to km/s
            float ax= -1.5*Omega*p->acc[i][0];
            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                    ax -= p->dx1[i][0]*q1 + p->dx2[i][0]*q2;
                    break;
                case FORCE_MODE_COLA1:
                    ax -= p->dx1[i][0]*q1;
                    break;
            } 
            po->v[i][d] = vfac*(p->v[i][d] + ax*dda);

            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                    po->v[i][d] += vfac * (p->dx1[i][d]*Dv + p->dx2[i][d]*Dv2);
                break;
                case FORCE_MODE_COLA1:
                    po->v[i][d] += vfac * (p->dx1[i][d]*Dv);
                break;
            }

            // Drift
            po->x[i][d] = p->x[i][d]; 
            if(FORCE_MODE & FORCE_MODE_PM) {
                po->x[i][d] += p->v[i][d]*dyyy;

            } 
            switch(FORCE_MODE) {
                case FORCE_MODE_COLA:
                case FORCE_MODE_2LPT:
                    po->x[i][d] += p->dx1[i][d]*da1 + p->dx2[i][d]*da2;
                break; 
                case FORCE_MODE_COLA1:
                case FORCE_MODE_ZA:
                    po->x[i][d] += p->dx1[i][d]*da1;
                break;
            }
        }

        po->id[i] = p->id[i];
    }

    po->np = np;
}

double growth_int(double a, void *param)
{
    return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


double growth(double a)
{
    int WORKSIZE = 100000;
    double hubble_a;

    hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda);

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;

    gsl_integration_qag(&F, 0, a, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return hubble_a * result;
}

double GrowthFactor(double astart, double aend)
{
    return growth(aend) / growth(astart);
}

