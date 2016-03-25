#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include <fastpm/cosmology.h>

static double growth_int(double a, void *param)
{
    double * p = (double*) param;
    double OmegaM = p[0];
    double OmegaLambda = p[1];
    return pow(a / (OmegaM + (1 - OmegaM - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


static double growth(double a, Cosmology c)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */

    int WORKSIZE = 100000;
    double hubble_a;

    hubble_a = Qfactor(a, c) / (a * a * a);

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {c.OmegaM, c.OmegaLambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return hubble_a * result;
}

/*
static double decayD(double a, Cosmology c){ // D_{-}, the decaying mode
    return sqrt(c.OmegaM/(a*a*a)+1.0-c.OmegaM);
}
*/
static double growthD(double a, Cosmology c) { // growth factor for LCDM
    return growth(a, c) / growth(1.0, c);
}

static double growthD2temp(double a, Cosmology c){
    double d = growthD(a, c);
    return d*d*pow(OmegaA(a, c), -1.0/143.);
}

double OmegaA(double a, Cosmology c) {
    return c.OmegaM/(c.OmegaM + (c.OmegaLambda)*a*a*a);
}
double DprimeQ(double a, double nGrowth, Cosmology c)
{
    /* This could have been Omega^(5/9) * Q / a * D1 for LCDM */ 
    // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
    double d = GrowthFactor(a, c);
    double r1 = Qfactor(a, c)*(d/a) * pow(OmegaA(a, c), 5.0 / 9);
/*
    double nDecay = 0.0;// not interested in decay modes in this code.
    double Nn = 1.0 / growth(1.0, c);
    double r2 = (  pow(decayD(a, c), -1.0 + nDecay)
            * pow(growthD(a, c),-1.0 + nGrowth)
            * (nGrowth*Nn - (3.0*(nDecay + nGrowth)*c.OmegaM *growthD(a, c))/(2.*a)));
//    printf("r1 = %g r2 = %g\n", r1, r2);
*/
    return r1;
}

double 
GrowthFactor(double a, Cosmology c)
{
    return growthD(a, c);
}

double GrowthFactor2(double a, Cosmology c) {// Second order growth factor
    return growthD2temp(a, c)/growthD2temp(1.0, c); // **???
}


double GrowthFactor2v(double a, Cosmology c){ // explanation is in main()
    /* This mess needs to be cleaned up. The extra pow is to match up
     * the original cola factor, since we no longer absorb D20 into dx2; 
     * D20 is absorbed to GrowthFactor2. */
    double d2= GrowthFactor2(a, c);
    return Qfactor(a, c)*(d2/a)*2.0
         * pow(OmegaA(a, c), 6.0/11.);
}

double Qfactor(double a, Cosmology c) {
 // Q\equiv a^3 H(a)/H0.
    return sqrt(c.OmegaM/(a*a*a)+c.OmegaLambda)*a*a*a;
}

#ifdef TEST_COSMOLOGY
int main() {
    /* the old COLA growthDtemp is 6 * pow(1 - c.OmegaM, 1.5) times growth */
    double a;
    Cosmology c = {
        .OmegaM = 0.3,
        .OmegaLambda = 0.7
    };

    for(c.OmegaM = 0.1; c.OmegaM < 0.6; c.OmegaM += 0.1) {
        double f = 6 * pow(1 - c.OmegaM, 1.5);
        c.OmegaLambda = 1 - c.OmegaM;
        double a = 0.1;
        printf("%g %g %g \n",
            c.OmegaM, 
                growth(a, c),
                GrowthFactor(a, c)
            );
    
    }
}

#endif
