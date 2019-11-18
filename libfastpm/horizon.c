#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>

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
        FastPMGrowthInfo gi;
        fastpm_growth_info_init(&gi, a, horizon->cosmology);
        horizon->growthfactor_a[i] = gi.D1;
    }
}

void
fastpm_horizon_destroy(FastPMHorizon * horizon)
{
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

double
HorizonGrowthFactor(double a, FastPMHorizon * horizon)
{
    double x = a * (horizon->size - 1);
    int l = floor(x);
    int r = l + 1;
    if(r >= horizon->size) {
        return horizon->growthfactor_a[horizon->size - 1];
    }
    if(l <= 0) {
        return horizon->growthfactor_a[0];
    }
    return horizon->growthfactor_a[l] * (r - x)
         + horizon->growthfactor_a[r] * (x - l);
}

void *
fastpm_horizon_solve_start()
{
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    return gsl_root_fsolver_alloc(T);
}

void
fastpm_horizon_solve_end(void * context)
{
    gsl_root_fsolver_free(context);
}

int
fastpm_horizon_solve(FastPMHorizon * horizon,
    void * context,
    double * solution,
    double a_i, double a_f,
    double (*func)(double a, void * userdata),
    void * userdata)
{

    int status;
    int iter = 0, max_iter;
    double r, x_lo=a_i, x_hi=a_f, eps;

    /* Reorganize to struct later */
    max_iter = 20;
    eps = 1e-5;

    gsl_function F;

    F.function = func;
    F.params = userdata;

    status = gsl_root_fsolver_set(context, &F, x_lo, x_hi);

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

        status = gsl_root_fsolver_iterate(context);
        r = gsl_root_fsolver_root(context);

        x_lo = gsl_root_fsolver_x_lower(context);
        x_hi = gsl_root_fsolver_x_upper(context);

        status = gsl_root_test_interval(x_lo, x_hi, eps, 0.0);
        //
        //Debug printout #2
        //fastpm_info("%5d [%.7f, %.7f] %.7f %.7f %.7f\n", iter, x_lo, x_hi, r, funct(r, &params), x_hi - x_lo);
        //

        if(status == GSL_SUCCESS || iter == max_iter ) {
            *solution = r;
            //
            // Debug printout #3.1
            //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 1);
            //
            return 1;
        }
    }
    while (status == GSL_CONTINUE);
    //
    // Debug printout #3.2
    //fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 0);
    //

    return 0;

}
