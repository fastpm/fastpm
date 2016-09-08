#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/cosmology.h>
#include <fastpm/lightcone.h>
#include <fastpm/logging.h>

void
fastpm_lc_init(FastPMLightCone * lc, Cosmology CP, size_t np_upper)
{
    /* Allocation */

    int size = 8192;

    lc->EventHorizonTable.size = size;
    lc->EventHorizonTable.Dc = malloc(sizeof(double) * size);

    int i;

    for(i = 0; i < lc->EventHorizonTable.size; i ++) {
        double a = 1.0 * i / (lc->EventHorizonTable.size - 1);
        lc->EventHorizonTable.Dc[i] = ComovingDistance(a, CP);
    }
}

void fastpm_lc_destroy(FastPMLightCone * lc)
{
    /* Free */

    free(lc->EventHorizonTable.Dc);
}

double fastpm_lc_horizon(FastPMLightCone * lc, double a)
{
    /* It may be worth to switch to log_a interpolation, but it only matters
     * at very high z (~ z = 9). */

    double x = a * (lc->EventHorizonTable.size - 1);
    int l = floor(x);
    int r = l + 1;
    if(r >= lc->EventHorizonTable.size) {
        return lc->EventHorizonTable.Dc[lc->EventHorizonTable.size - 1];
    }
    if(l <= 0) {
        return lc->EventHorizonTable.Dc[0];
    }
    return lc->EventHorizonTable.Dc[l] * (r - x)
         + lc->EventHorizonTable.Dc[r] * (x - l);
}

double funct(double x, void *params)
{
  struct funct_params *Fp = (struct funct_params *) params;

  double a = Fp->a;
  double b = Fp->b;
  FastPMLightCone *lc = Fp->lc;

  return (a * x + b) - fastpm_lc_horizon(lc, x);
}

//int fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p, int i, double * solution)
int fastpm_lc_intersect(FastPMLightCone * lc, double * solution, double a, double b)
{
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	int status;
	int iter = 0, max_iter;
	double r, x_lo, x_hi, eps;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);

	gsl_function F;
	
	/* Reorganize to struct later */
	x_lo = 0.0;
	x_hi = 1.0;
	max_iter = 100;
	eps = 0.0000001;
	struct funct_params params = {lc, a, b};
	
	F.function = &funct;
	F.params = &params;
	
	gsl_set_error_handler_off(); // Turn off GSL error handler
	status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);
	
	if(status == GSL_EINVAL || status == GSL_EDOM) { 
		// Debug printout #0
		fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 0);
		//
		return 0; 
	} // Error in value or out of range
	
	do
	{
		iter++;
		
		// Debug printout #1
		if(iter == 1) {
		fastpm_info("ID | [x_lo, x_hi] | r | funct(r) | x_hi - x_lo\n");
		}
		//
		
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		
		status = gsl_root_test_interval(x_lo, x_hi, eps, 0.0);
		
		//Debug printout #2
		fastpm_info("%5d [%.7f, %.7f] %.7f %.7f %.7f\n", iter, x_lo, x_hi, r, funct(r, &params), x_hi - x_lo);
		//
		
		if(status == GSL_SUCCESS) {
			solution = &r;
			gsl_root_fsolver_free(s);
			
			// Debug printout #3.1
			fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 1);
			//
			return 1;
		}
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	// Debug printout #3.2
	fastpm_info("fastpm_lc_intersect() called with parameters %.7f and %.7f, returned status %d.\n\n", a, b, 0);
	//

    return 0;
}

#if 0

void fastpm_lc_write(FastPMLightCone * lc, const char * filename)
{

}

double _unknown(double a, void * params)
{

    double xo[3];
    double horizon;

    horizon = fastpm_lc_horizon(lc, a);

    fastpm_drift_one(drift, p, i, xo, a);

    return horizon - xo[2];
}

int fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p, int i, double * solution)
{

    /* XXX: culling */

    gsl_fsolve(solution);

    if( *solution >= drift->a_i
     && *solution <= drift->a_f)
    {
        return 1;
    }
    return 0;
}

int fastpm_lc_scan(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p)
{
    for(i in p) {
        status = fastpm_lc_intersect(lc, drift, p, i, solution)
        if(dint is good)
            add(i, lc);
    }
}
#endif
