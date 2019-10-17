#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

//include Fermi-Dirac integration table for neutrinos
#include <fastpm/Ftable.h>

void fastpm_fd_interp_init(FastPMFDInterp * FDinterp)
{
    FDinterp->size = Fsize;  //might want to neaten how Fsize comes into play 
    FDinterp->F   = gsl_interp_alloc(gsl_interp_linear, Fsize);
    FDinterp->DF  = gsl_interp_alloc(gsl_interp_linear, Fsize);
    FDinterp->DDF = gsl_interp_alloc(gsl_interp_linear, Fsize);
    FDinterp->acc = gsl_interp_accel_alloc();
    
    //FIXME: do i want to init here? I think so.
    gsl_interp_init(FDinterp->F,   Ftable[0], Ftable[1], Fsize);
    gsl_interp_init(FDinterp->DF,  Ftable[0], Ftable[2], Fsize);
    gsl_interp_init(FDinterp->DDF, Ftable[0], Ftable[3], Fsize);
}

double fastpm_do_fd_interp(FastPMFDInterp * FDinterp, int F_id, double y)
{
    /* Gets the interpolated value of Ftable[F_id] at y
       F_id: 1 for F, 2 for F', 3 for F'' */
    double res = 0;  //give this a value to make compiler happy?
    switch (F_id){
        case 1:
            res = gsl_interp_eval(FDinterp->F, Ftable[0], Ftable[F_id], y, FDinterp->acc);
        break;
        case 2:
            res = gsl_interp_eval(FDinterp->DF, Ftable[0], Ftable[F_id], y, FDinterp->acc);
        break;
        case 3:
            res = gsl_interp_eval(FDinterp->DDF, Ftable[0], Ftable[F_id], y, FDinterp->acc);
        break;
        default:
            fastpm_raise(-1, "Wrong id for FD table.\n");
    }
    gsl_interp_accel_reset(FDinterp->acc);   // reset accelerator for next time
    return res;
}

void fastpm_fd_interp_free(FastPMFDInterp * FDinterp)
{
    gsl_interp_free(FDinterp->F);
    gsl_interp_free(FDinterp->DF);
    gsl_interp_free(FDinterp->DDF);
    gsl_interp_accel_free(FDinterp->acc);
}
