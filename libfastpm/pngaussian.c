#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"

static double
fastpm_png_potential(double k, FastPMPNGaussian * png)
{
    /*FIXME: use the correct primordial power to construct the potential !*/
    return 1.0 / sqrt(png->Volume);
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
    for(i = 0; i < pm_size(pm); i ++) {
        g_x[i] = g_x[i] + png->fNL * g_x[i] * g_x[i];
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
