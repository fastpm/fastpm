#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>

double HubbleDistance = 2997.92458; /* Mpc/h */
double HubbleConstant = 100.0; /* Mpc/h / km/s*/

static double growth_int(double a, void *param)
{
    double * p = (double*) param;
    double OmegaM = p[0];
    double OmegaLambda = p[1];
    return pow(a / (OmegaM + (1 - OmegaM - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


static double growth(double a, FastPMCosmology * c)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {c->OmegaM, c->OmegaLambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return HubbleEa(a, c) * result;
}

double OmegaA(double a, FastPMCosmology * c) {
    return c->OmegaM/(c->OmegaM + (c->OmegaLambda)*a*a*a);
}

double GrowthFactor(double a, FastPMCosmology * c) { // growth factor for LCDM
    return growth(a, c) / growth(1.0, c);
}

double DLogGrowthFactor(double a, FastPMCosmology * c) {
    /* Or OmegaA^(5/9) */
    return pow(OmegaA(a, c), 5.0 / 9);
}

double GrowthFactor2(double a, FastPMCosmology * c) {
    /* Second order growth factor */
    /* 7 / 3. is absorbed into dx2 */
    double d = GrowthFactor(a, c);
    return d * d * pow(OmegaA(a, c) / OmegaA(1.0, c), -1.0/143.);
}

double DLogGrowthFactor2(double a, FastPMCosmology * c) {
    return 2 * pow(OmegaA(a, c), 6.0/11.);
}

double HubbleEa(double a, FastPMCosmology * c)
{
    /* H(a) / H0 */
    return sqrt(c->OmegaM/(a*a*a)+c->OmegaLambda);
}
double DHubbleEaDa(double a, FastPMCosmology * c) {
    /* d E / d a*/
    double E = HubbleEa(a, c);
    return 0.5 / E * (-3 * c->OmegaM / (a * a * a * a));
}
double D2HubbleEaDa2(double a, FastPMCosmology * c) {
    double E = HubbleEa(a, c);
    double dEda = DHubbleEaDa(a, c);
    return - dEda * dEda / E + dEda * (-4 / a);
}
double DGrowthFactorDa(double a, FastPMCosmology * c) {
    double E = HubbleEa(a, c);

    double EI = growth(1.0, c);

    double t1 = DHubbleEaDa(a, c) * GrowthFactor(a, c) / E;
    double t2 = E * pow(a * E, -3) / EI;
    return t1 + t2;
}
double D2GrowthFactorDa2(double a, FastPMCosmology * c) {
    double d2Eda2 = D2HubbleEaDa2(a, c);
    double dEda = DHubbleEaDa(a, c);
    double E = HubbleEa(a, c);
    double EI = growth(1.0, c);
    double t1 = d2Eda2 * GrowthFactor(a, c) / E;
    double t2 = (dEda + 3 / a * E) * pow(a * E, -3) / EI;
    return t1 - t2;
}

static double
comoving_distance_int(double a, void * params)
{
    FastPMCosmology * c = (FastPMCosmology * ) params;
    return 1 / (a * a * HubbleEa(a, c));
}

/* In Hubble Distance */
double ComovingDistance(double a, FastPMCosmology * c) {

    /* We tested using ln_a doesn't seem to improve accuracy */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &comoving_distance_int;
    F.params = (void*) c;

    gsl_integration_qag(&F, a, 1, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41,
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return result;
}

void
fastpm_horizon_init(FastPMHorizon * horizon, FastPMCosmology * cosmology)
{
    gsl_set_error_handler_off(); // Turn off GSL error handler

    horizon->cosmology = cosmology;
    horizon->size = 8192;
    horizon->da = 1.0 / (horizon->size - 1);
    int i;
    for (i = 0; i < horizon->size; i ++) {
        double a = 1.0 * i / (horizon->size - 1);
        horizon->xi_a[i] = HubbleDistance * ComovingDistance(a, horizon->cosmology);
    }
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    horizon->gsl = gsl_root_fsolver_alloc(T);
}

void
fastpm_horizon_destroy(FastPMHorizon * horizon)
{
    gsl_root_fsolver_free(horizon->gsl);
}

double
HorizonDistance(double a, FastPMHorizon * horizon)
{
    double x = a * (horizon->size - 1);
    int l = floor(x);
    int r = l + 1;
    if(r >= horizon->size) {
        return horizon->xi_a[horizon->size - 1];
    }
    if(l <= 0) {
        return horizon->xi_a[0];
    }
    return horizon->xi_a[l] * (r - x)
         + horizon->xi_a[r] * (x - l);
}

int
fastpm_horizon_solve(FastPMHorizon * horizon,
    double * solution,
    double a_i, double a_f,
    double (*func)(double a, void * userdata),
    void * userdata)
{

    int status;
    int iter = 0, max_iter;
    double r, x_lo=a_i, x_hi=a_f, eps;

    /* Reorganize to struct later */
    max_iter = 100;
    eps = 1e-7;

    gsl_function F;

    F.function = func;
    F.params = userdata;

    status = gsl_root_fsolver_set(horizon->gsl, &F, x_lo, x_hi);

    if(status == GSL_EINVAL || status == GSL_EDOM) {
        /** Error in value or out of range **/
        return 0;
    }

    do
    {
        iter++;
        //
        // Debug printout #1
        //if(iter == 1) {
        //fastpm_info("ID | [x_lo, x_hi] | r | funct(r) | x_hi - x_lo\n");
        //}
        //

        status = gsl_root_fsolver_iterate(horizon->gsl);
        r = gsl_root_fsolver_root(horizon->gsl);

        x_lo = gsl_root_fsolver_x_lower(horizon->gsl);
        x_hi = gsl_root_fsolver_x_upper(horizon->gsl);

        status = gsl_root_test_interval(x_lo, x_hi, eps, 0.0);
        //
        //Debug printout #2
        //fastpm_info("%5d [%.7f, %.7f] %.7f %.7f %.7f\n", iter, x_lo, x_hi, r, funct(r, &params), x_hi - x_lo);
        //

        if(status == GSL_SUCCESS) {
            *solution = r;
            //
            // Debug printout #3.1
            //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 1);
            //
            return 1;
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    //
    // Debug printout #3.2
    //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 0);
    //

    return 0;

}

#ifdef TEST_COSMOLOGY
int main() {
    /* the old COLA growthDtemp is 6 * pow(1 - c.OmegaM, 1.5) times growth */
    double a;
    FastPMCosmology c[1] = {{
        .OmegaM = 0.3,
        .OmegaLambda = 0.7
    }};

    printf("OmegaM D dD/da d2D/da2 D2 E dE/dA d2E/da2 \n");
    for(c->OmegaM = 0.1; c->OmegaM < 0.6; c->OmegaM += 0.1) {
        double f = 6 * pow(1 - c->OmegaM, 1.5);
        c->OmegaLambda = 1 - c->OmegaM;
        double a = 0.8;
        printf("%g %g %g %g %g %g %g %g %g\n",
            c->OmegaM, 
            ComovingDistance(a, c),
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
