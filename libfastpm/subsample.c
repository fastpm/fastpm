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

/* Computes particle volume number density [1 / (Mpc/h)^3] 
 * to reach the ell_lim resolution at given redshift */
double VolumeDensityFromEll(double ell_lim, double z, FastPMCosmology * c)
{
    double theta_lim = M_PI / ell_lim;
    double scale_fac = 1 / (1 + z);
    double r = ComovingDistance(scale_fac,c) * HubbleDistance;
    double s_lim = r * theta_lim;
    double res_lim = pow(1 / s_lim, 3);
    return res_lim;
}

