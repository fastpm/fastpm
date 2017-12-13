#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"

static double
fastpm_png_potential(double k, FastPMPNGaussian * png)
{
    /* Returns primordial power spectrum P_Phi(k), assuming k is in
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

    /*
     *  The parameterization is such that
     *
     *  <d**2, d> = 4 * sigma_L ** 2 * fNL <d, phi>
     *            = 4 * sigma_l ** 2 * fNL * <d, d> / M
     *   
     *  where M = 2 / 3 (D_H k) ** 2 T(k) D_LSS(z) / Om
     *
     *  where T(k) -> 1 when k -> 0
     *
     *  where D(z)(1+z) = 1 in matter era and D_H is c / H_0, hubble distance.
     *
     * */
    if (k == 0) return 0.0;

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

    double pk = png->pkfunc(k, png->pkdata);

    /* powerspec = transfer^2 * pot, so we remove pot */
    double transfer = sqrt(pk / fastpm_png_potential(k, png));
    return transfer;
}

static void
fastpm_png_transform_potential(PM * pm, FastPMFloat * delta_k, FastPMPNGaussian * png)
{
    FastPMFloat * g_x = pm_alloc(pm);
    FastPMFloat * g_x2 = pm_alloc(pm);

    /* The first order, gaussian piece is the full spetram*/
    pm_assign(pm, delta_k, g_x);
    pm_c2r(pm, g_x);

    /* The second order, NG piece must be truncated to avoid Dirac folding. */
    pm_assign(pm, delta_k, g_x2);
    /* MS: Zero-pad/truncate high k to avoid spurious Dirac delta images. */
    fastpm_apply_lowpass_transfer(pm, g_x2, g_x2, png->kmax_primordial);
    pm_c2r(pm, g_x2);

    ptrdiff_t i;
    double avg_g_squared = 0.0;
    PMXIter xiter;

    for(pm_xiter_init(pm, &xiter);
       !pm_xiter_stop(&xiter);
        pm_xiter_next(&xiter)) {

        avg_g_squared += g_x2[xiter.ind] * g_x2[xiter.ind];
    }

    MPI_Allreduce(MPI_IN_PLACE, &avg_g_squared, 1, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    avg_g_squared /= pm_norm(pm);

    /* Expect <Phi**2(x)> = \int dq/(2 pi**2) q**2 P_Phi(q) */

    /* FIXME:
     * This is off because the large scale modes with large power are sampled poorly.
     * in DFT; see comment above for testing the normalization factors. */
    double avg_g_squared_exp = 0.0;
    double q, dq;
    dq = png->kmax_primordial/1e6;
    double k0 = 2 * M_PI / pm->BoxSize[0];
    for(q = k0; q < png->kmax_primordial ; q += dq) {
        double p = fastpm_png_potential(q, png);
        avg_g_squared_exp = avg_g_squared_exp + q*q * p;
    }
    avg_g_squared_exp *= dq/(2.0*M_PI*M_PI);

    fastpm_info("Expected_avg_g_squared: %g, when there is no gaussian variance\n", avg_g_squared_exp);
    fastpm_info("avg_g_squared: %g, %g\n", avg_g_squared, avg_g_squared*avg_g_squared);

    /* Use the realization mean because we want the mean of linear field to be exactly zero. */
    for(i = 0; i < pm_allocsize(pm); i ++) {
        g_x[i] = g_x[i] + png->fNL * ( g_x2[i] * g_x2[i] - avg_g_squared );
    }

    /* MS: Print some info */
    fastpm_info("Induced PNG with fNL=%g g_x[0] = %g\n", png->fNL, g_x[0]);

    pm_r2c(pm, g_x, delta_k);

    pm_free(pm, g_x2);
    pm_free(pm, g_x);
}

void
fastpm_png_induce_correlation(FastPMPNGaussian * png, PM * pm, FastPMFloat * delta_k)
{
    png->Volume = pm->Volume;

    fastpm_ic_induce_correlation(pm, delta_k, (fastpm_fkfunc) fastpm_png_potential, png);

    fastpm_png_transform_potential(pm, delta_k, png);

    fastpm_apply_any_transfer(pm, delta_k, delta_k, (fastpm_fkfunc) fastpm_png_transfer_function, png);
}

/* vim: set ts=4 sw=4 sts=4 expandtab */
