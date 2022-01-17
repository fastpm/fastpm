//
//  subsampling.c
//  
//
//  Created by yici zhong on 2021/12/23.
//

//#astropy.cosmology in c library: FlatLambdaCDM,comoving_volume, etc.

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include <fastpm/libfastpm.h>

#define vlight 3.e5

double FlatLambdaCDM_kpc_per_arcmin(double l)
{
    return l / 3.4377468;
}

static
double degrees(double x)
{
    return x*180/M_PI;
}

//function to calculate subsampling rate
double calc_subrate(double ell_lim, double z, double res_box, FastPMCosmology * c)
{
    double theta_lim_rad = M_PI / ell_lim;
    double theta_lim_arcmin = degrees(theta_lim_rad) * 60;
    //printf("%f, arcmin(basic)\n",theta_lim_arcmin);
    double scale_fac = 1 / (1 + z);
    double com_dis = ComovingDistance(scale_fac,c)*vlight/(c->h*100);
    double kpc_per_arcmin = FlatLambdaCDM_kpc_per_arcmin(com_dis);
    double Mpc_res_lim = kpc_per_arcmin * theta_lim_arcmin / 1000.;
    double res_lim = pow(1 / Mpc_res_lim, 3);
    double rate_sub = res_lim / res_box;
    /* FIXME: in principle we can replicate particles to achieve a rate > 1.
     * probably want to move this clipping to the caller side.*/
    if (rate_sub > 1) {
        rate_sub = 1; //for low redshift, no subsample
    }
    //fprintf(f2, "%7.6lf, %10.7lf\n",z, rate_sub);
    return rate_sub;
}

#undef vlight
