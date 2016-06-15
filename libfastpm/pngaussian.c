#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
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
     *    scalar_pivot: 
     *      Pivot scale k_pivot in 1/Mpc. Same as pivot_scalar in CAMB.
     *      Example: 0.05 [1/Mpc].
     *
     *    scalar_spectral_index: 
     *      Primordial spectral index n_s. Same as scalar_spectral_index
     *      in CAMB. Example: 0.9653.
     */

    /* MS: More notes on primordial non-Gaussianity parameters
     *
     *   kmax_primordial_over_knyquist:
     *     To avoid spurious Dirac foldings when computing Phi^2(x) on a grid, 
     *     truncate (zero-pad) the primordial Phi power spectrum at 
     *     k>=kmax_primordial, where
     *    
     *       kmax_primordial = kmax_primordial_over_knyquist * knyquist
     *
     *     and knyquist = N_grid/2 * 2 pi/boxsize. Example: 0.25 will set
     *     Phi(k)=0 for k>=knyquist/4. Should be less than 0.5 or 0.25.
     */

    if (k == 0) return 0.0;

    /* MS: Zero-pad/truncate high k to avoid spurious Dirac delta images. */
    if (k >= png->kmax_primordial) return 0.0;

    double k_pivot_in_h_over_Mpc, P_Phi_k;

    /* k_pivot is in 1/Mpc, so need to divide by h to get it in h/Mpc.  */
    k_pivot_in_h_over_Mpc = png->scalar_pivot / png->h;

    /* Compute A_s / k^3 */
    P_Phi_k = png->scalar_amp;

    /* disable slope and tilt to use a flat transfer function testing numerical normalization factors */

    /* Slope */
    P_Phi_k /= k*k*k;

    /* Tilt */
    P_Phi_k *= pow(k/k_pivot_in_h_over_Mpc, png->scalar_spectral_index - 1.0);

    /* Prefactor */
    P_Phi_k *= 9.0/25.0 * 2.0 * M_PI * M_PI;


    return P_Phi_k;
}

static double
fastpm_png_transfer_function(double k, FastPMPNGaussian * png)
{
    if (k == 0) return 0.0;

    /* MS: Zero-pad/truncate high k to avoid spurious Dirac delta images. */
    if (k >= png->kmax_primordial) return 0.0;

    double transfer = sqrt(png->pkfunc(k, png->pkdata));

    /* powerspec = transfer^2 * pot, so we remove pot */
    transfer /= fastpm_png_potential(k, png);

    return transfer;
}

static void
fastpm_png_transform_potential(PM * pm, FastPMFloat * g_x, FastPMPNGaussian * png)
{
    ptrdiff_t i;
    double avg_g_squared = 0.0;
    /* MS: Should we better do this globally over all processors? */
    PMXIter xiter;

    for(pm_xiter_init(pm, &xiter);
       !pm_xiter_stop(&xiter);
        pm_xiter_next(&xiter)) {

        avg_g_squared += g_x[xiter.ind] * g_x[xiter.ind];

    }

    MPI_Allreduce(MPI_IN_PLACE, &avg_g_squared, 1, MPI_DOUBLE, MPI_SUM, pm->Comm2D);

    avg_g_squared /= pm_norm(pm);

    for(i = 0; i < pm_size(pm); i ++) {
        g_x[i] = g_x[i] + png->fNL * ( g_x[i] * g_x[i] - avg_g_squared );
    }

    /* MS: Print some info */
    fastpm_info("Induced PNG with fNL=%g\n",png->fNL);
    fastpm_info("avg_g_squared: %g, %g\n", avg_g_squared, avg_g_squared*avg_g_squared);

    /* Expect <Phi**2(x)> = \int dq/(2 pi**2) q**2 P_Phi(q) */

    /* FIXME:
     * This is off because the large scale modes with large power are sampled poorly.
     * in DFT; see comment above for testing the normalization factors. */
    double avg_g_squared_exp = 0.0;
    double q, dq;
    dq = png->kmax_primordial/1e6;
    double k0 = 2 * M_PI / pm->BoxSize[0];
    for(q = k0; q < png->kmax_primordial ; q += dq) {
        double t = fastpm_png_potential(q,png);
        avg_g_squared_exp = avg_g_squared_exp + q*q *(t * t);
    }
    avg_g_squared_exp *= dq/(2.0*M_PI*M_PI);
    fastpm_info("Expected_avg_g_squared: %g, assuming no gaussian variance\n", avg_g_squared_exp);

}

void
fastpm_png_induce_correlation(FastPMPNGaussian * png, PM * pm, FastPMFloat * delta_k)
{
    FastPMFloat * g_x = pm_alloc(pm);
    png->Volume = pm->Volume;

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) fastpm_png_potential, png);
    fastpm_apply_multiply_transfer(pm, delta_k, delta_k, 1 / sqrt(png->Volume));

    pm_assign(pm, delta_k, g_x);
    pm_c2r(pm, g_x);

    fastpm_png_transform_potential(pm, g_x, png);

    pm_r2c(pm, g_x, delta_k);
    pm_free(pm, g_x);

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) fastpm_png_transfer_function, png);
    fastpm_apply_multiply_transfer(pm, delta_k, delta_k, 1 / sqrt(png->Volume));
}

/* vim: set ts=4 sw=4 sts=4 expandtab */
