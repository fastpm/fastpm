#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"

static double
fastpm_png_potential(double k, FastPMPNGaussian * png)
{
    /* Returns primordial power spectrum P_Phi(k)/volume, assuming k is in
     * h/Mpc units.
     */

    /* MS: Initial power spectrum notes.
     *
     *  Use the same standard running power law model for the primordial super- 
     *  horizon power spectrum of curvature perturbations as CAMB, i.e. 
     *  
     *    k^3 P_zeta(k) / (2 pi^2) = A_s (k/k_pivot)^(n_s-1), 
     *  
     *  where A_s is the amplitude at the pivot scale k=k_pivot. The power of the 
     *  primordial potential Phi is therefore 
     *   
     *    P_Phi(k) = (9/25) P_zeta(k) 
     *             = (9/25) (2 pi^2 / k^3) A_s (k/k_pivot)^(n_s-1). 
     *
     *  In parameter file, A_s, k_pivot and n_s are specified as follows:
     *
     *    scalar_amp: 
     *      Amplitude A_s of primordial power spectrum at pivot scale.
     *      Same as scalar_amp in CAMB. Example: 2.130624e-9.
     *
     *    k_pivot: 
     *      Pivot scale k_pivot in 1/Mpc. Same as pivot_scalar in CAMB.
     *      Example: 0.05 [1/Mpc].
     *
     *    scalar_spectral_index: 
     *      Primordial spectral index n_s. Same as scalar_spectral_index
     *      in CAMB. Example: 0.9653.
     */

    /*FIXME: use the correct primordial power to construct the potential !*/    

    double k_pivot_in_h_over_Mpc, P_Phi_k;

    /* k_pivot is in 1/Mpc, so need to divide by h to get it in h/Mpc.  */
    k_pivot_in_Mpc_over_h = png->k_pivot / png->h;

    /* Compute A_s / k^3 */
    P_Phi_k = png->scalar_amp / k;
    P_Phi_k /= k*k;
    /* Prefactor */
    P_Phi_k *= 9.0/25.0 * 2.0 * M_PI * M_PI;
    /* Tilt */
    P_Phi_k *= pow(k/k_pivot_in_Mpc_over_h, png->scalar_spectral_index - 1.0)

    return P_Phi_k / sqrt(png->Volume);
}

static double
fastpm_png_transfer_function(double k, FastPMPNGaussian * png)
{
    double transfer = sqrt(png->pkfunc(k, png->pkdata));
    /* powerspec = transfer^2 * pot, so we remove pot */
    transfer /= fastpm_png_potential(k, png);
    /* don't forget the volume factor */
    transfer *= 1.0 / sqrt(png->Volume);
    return transfer;
}

static void
fastpm_png_transform_potential(PM * pm, FastPMFloat * g_x, FastPMPNGaussian * png)
{
    /*FIXME: transform the potential !*/
    ptrdiff_t i;
    double avg_g_squared = 0.0;
    for(i = 0; i < pm_size(pm); i ++) {
	avg_g_squared += g_x[i] * g_x[i];
    }
    for(i = 0; i < pm_size(pm); i ++) {
        g_x[i] = g_x[i] + png->fNL * ( g_x[i] * g_x[i] - avg_g_squared );
    }
}

void
fastpm_png_induce_correlation(FastPMPNGaussian * png, PM * pm, FastPMFloat * delta_k)
{
    FastPMFloat * g_x = pm_alloc(pm);
    png->Volume = pm->Volume;

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) fastpm_png_potential, png);

    pm_assign(pm, delta_k, g_x);
    pm_c2r(pm, g_x);

    fastpm_png_transform_potential(pm, g_x, png);

    pm_r2c(pm, g_x, delta_k);
    pm_free(pm, g_x);

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) fastpm_png_transfer_function, png);
}
