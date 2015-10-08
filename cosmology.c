#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

static double growth_int(double a, void *param)
{
    double * p = (double*) param;
    double Omega = p[0];
    double OmegaLambda = p[1];
    return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


static double growth(double a, double Omega, double OmegaLambda)
{
    int WORKSIZE = 100000;
    double hubble_a;

    hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda);

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {Omega, OmegaLambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return hubble_a * result;
}

double 
GrowthFactor(double astart, double aend, double Omega, double OmegaLambda)
{
    return growth(aend, Omega, OmegaLambda) / growth(astart, Omega, OmegaLambda);
}

