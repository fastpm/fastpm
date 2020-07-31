#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

//include Fermi-Dirac integration table for ncdm
#include <fastpm/Ftable.h>

void fastpm_fd_interp_init(FastPMFDInterp * FDinterp)
{
    FDinterp->size = Fsize;  // FIXME: might want to neaten how Fsize comes into play 
    FDinterp->F   = gsl_interp_alloc(gsl_interp_cspline, Fsize);
    FDinterp->DF  = gsl_interp_alloc(gsl_interp_cspline, Fsize);
    FDinterp->DDF = gsl_interp_alloc(gsl_interp_cspline, Fsize);
    FDinterp->acc = gsl_interp_accel_alloc();
    
    gsl_interp_init(FDinterp->F,   Ftable[0], Ftable[1], Fsize);
    gsl_interp_init(FDinterp->DF,  Ftable[0], Ftable[2], Fsize);
    gsl_interp_init(FDinterp->DDF, Ftable[0], Ftable[3], Fsize);
}

double fastpm_do_fd_interp(FastPMFDInterp * FDinterp, int F_id, double y)
{
    /* Gets the interpolated value of Ftable[F_id] at y
       F_id: 1 for F, 2 for F', 3 for F'' */
    int status = 1;
    double res;
    switch (F_id){
        case 1:
            status = gsl_interp_eval_e(FDinterp->F, Ftable[0], Ftable[F_id], y, FDinterp->acc, &res);
        break;
        case 2:
            status = gsl_interp_eval_e(FDinterp->DF, Ftable[0], Ftable[F_id], y, FDinterp->acc, &res);
        break;
        case 3:
            status = gsl_interp_eval_e(FDinterp->DDF, Ftable[0], Ftable[F_id], y, FDinterp->acc, &res);
        break;
        default:
            fastpm_raise(-1, "Wrong id for FD table.\n");
    }
    if (status) {
        fastpm_raise(-1, "FD interpolation failed. Value outside of range.\n");
    }
    gsl_interp_accel_reset(FDinterp->acc);   // reset accelerator for next time
    return res;
}

void fastpm_fd_interp_destroy(FastPMFDInterp * FDinterp)
{
    gsl_interp_free(FDinterp->F);
    gsl_interp_free(FDinterp->DF);
    gsl_interp_free(FDinterp->DDF);
    gsl_interp_accel_free(FDinterp->acc);
}
