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

//function to calculate subsampling rate
double calc_subrate(double ell_lim, double z, double res_box, FastPMCosmology * c)
{
    double theta_lim = M_PI / ell_lim;
    double scale_fac = 1 / (1 + z);
    double r = ComovingDistance(scale_fac,c) * HubbleDistance;
    double s_lim = r * theta_lim;
    double res_lim = pow(1 / s_lim, 3);
    double rate_sub = res_lim / res_box;
    /* FIXME: in principle we can replicate particles to achieve a rate > 1.
     * probably want to move this clipping to the caller side.*/
    if (rate_sub > 1) {
        rate_sub = 1; // for low redshift, no subsample
    }
    return rate_sub;
}

