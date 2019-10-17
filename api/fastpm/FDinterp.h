#include <gsl/gsl_spline.h> //should this be in an h file?

typedef struct FastPMFDInterp{
    /* maybe would be nicer to index these, ratehr than name
      then you could loop thru nicer*/
    size_t size;
    gsl_interp * F;
    gsl_interp * DF;
    gsl_interp * DDF;
    gsl_interp_accel * acc;
    
} FastPMFDInterp;

void fastpm_fd_interp_init(FastPMFDInterp * FDinterp);

double fastpm_do_fd_interp(FastPMFDInterp * FDinterp, int F_id, double y);

void fastpm_fd_interp_free(FastPMFDInterp * FDinterp);
