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

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {c.OmegaM, c.OmegaLambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return HubbleEa(a, c) * result;
}

double OmegaA(double a, Cosmology c) {
    return c.OmegaM/(c.OmegaM + (c.OmegaLambda)*a*a*a);
}

double GrowthFactor(double a, Cosmology c) { // growth factor for LCDM
    return growth(a, c) / growth(1.0, c);
}

double DLogGrowthFactor(double a, Cosmology c) {
    /* Or OmegaA^(5/9) */
    return pow(OmegaA(a, c), 5.0 / 9);
}

double GrowthFactor2(double a, Cosmology c) {
    /* Second order growth factor */
    /* 7 / 3. is absorbed into dx2 */
    double d = GrowthFactor(a, c);
    return d * d * pow(OmegaA(a, c) / OmegaA(1.0, c), -1.0/143.);
}

double DLogGrowthFactor2(double a, Cosmology c) {
    return 2 * pow(OmegaA(a, c), 6.0/11.);
}

double HubbleEa(double a, Cosmology c)
{
    /* H(a) / H0 */
    return sqrt(c.OmegaM/(a*a*a)+c.OmegaLambda);
}
double DHubbleEaDa(double a, Cosmology c) {
    /* d E / d a*/
    double E = HubbleEa(a, c);
    return 0.5 / E * (-3 * c.OmegaM / (a * a * a * a));
}
double D2HubbleEaDa2(double a, Cosmology c) {
    double E = HubbleEa(a, c);
    double dEda = DHubbleEaDa(a, c);
    return - dEda * dEda / E + dEda * (-4 / a);
}
double DGrowthFactorDa(double a, Cosmology c) {
    double E = HubbleEa(a, c);

    double EI = growth(1.0, c);

    double t1 = DHubbleEaDa(a, c) * GrowthFactor(a, c) / E;
    double t2 = E * pow(a * E, -3) / EI;
    return t1 + t2;
}
double D2GrowthFactorDa2(double a, Cosmology c) {
    double d2Eda2 = D2HubbleEaDa2(a, c);
    double dEda = DHubbleEaDa(a, c);
    double E = HubbleEa(a, c);
    double EI = growth(1.0, c);
    double t1 = d2Eda2 * GrowthFactor(a, c) / E;
    double t2 = (dEda + 3 / a * E) * pow(a * E, -3) / EI;
    return t1 - t2;
}

#ifdef TEST_COSMOLOGY
int main() {
    /* the old COLA growthDtemp is 6 * pow(1 - c.OmegaM, 1.5) times growth */
    double a;
    Cosmology c = {
        .OmegaM = 0.3,
        .OmegaLambda = 0.7
    };

    printf("OmegaM D dD/da d2D/da2 D2 E dE/dA d2E/da2 \n");
    for(c.OmegaM = 0.1; c.OmegaM < 0.6; c.OmegaM += 0.1) {
        double f = 6 * pow(1 - c.OmegaM, 1.5);
        c.OmegaLambda = 1 - c.OmegaM;
        double a = 0.8;
        printf("%g %g %g %g %g %g %g %g\n",
            c.OmegaM, 
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
