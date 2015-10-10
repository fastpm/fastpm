/* libfastpm: */
#include <mpi.h>
#include <string.h>
#include "libfastpm.h"
#include "msg.h"

static int 
fastpm_particle_to_mesh(PM * pm, PMStore * p);

int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm) 
{
    msg_init(comm);
    msg_set_loglevel(verbose);

    pm_store_init(p);

    pm_store_alloc_evenly(p, pow(nc, 3), 2.0, comm);

    return 0;
}

static double 
tk_eh(double k, struct fastpm_powerspec_eh_params * params)		/* from Martin White */
{
    double q, theta, ommh2, a, s, gamma, L0, C0;
    double tmp;
    double omegam, ombh2, hubble;

    /* other input parameters */
    hubble = params->hubble_param;
    omegam = params->omegam;
    ombh2 = params->omegab * hubble * hubble;

    theta = 2.728 / 2.7;
    ommh2 = omegam * hubble * hubble;
    s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
    a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
        + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
    gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
    gamma *= omegam * hubble;
    q = k * theta * theta / gamma;
    L0 = log(2. * exp(1.) + 1.8 * q);
    C0 = 14.2 + 731. / (1. + 62.5 * q);
    tmp = L0 / (L0 + C0 * q * q);
    return (tmp);
}

double 
fastpm_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param)	/* Eisenstein & Hu */
{
    return param->Norm * k * pow(tk_eh(k, param), 2);
}

void 
fastpm_fill_deltak(PM * pm, real_t * deltak, int seed, fastpm_pkfunc pk, void * pkdata) 
{
    pm_start(pm);
    pm_ic_fill_gaussian_gadget(pm, seed, pk, pkdata);
    memcpy(deltak, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);
    pm_stop(pm);
}

void 
fastpm_evolve_2lpt(PM * pm, PMStore * pdata, 
        double a, double omega_m, 
        real_t * deltak_0, real_t * deltak_1, MPI_Comm comm) 
{

    pm_start(pm);

    memcpy(pm->canvas, deltak_0, sizeof(pm->canvas[0]) * pm->allocsize);

    pm_2lpt_main(pm, pdata, comm);

    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */
    int d;
    double shift[3];
    for(d = 0; d < 3; d ++){
        shift[d] = 0.0; //pm->CellSize[d] * 0.5;
    }

    /* predict particle positions by 2lpt */
    pm_2lpt_set_initial(a, pdata, shift, omega_m);

    /* paint to mesh */
    fastpm_particle_to_mesh(pm, pdata);

    pm_r2c(pm);

    /* copy out the results */
    memcpy(deltak_1, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);

    pm_stop(pm);
}

void 
fastpm_evolve_pm(PM * pm, PMStore * pdata, 
        double ainit, double afinal, int nsteps, double omega_m, 
        real_t * deltak_0, real_t * deltak_1, MPI_Comm comm) 
{
    fastpm_evolve_2lpt(pm, pdata, ainit, omega_m, deltak_0, deltak_1, comm);

    /* now do the steps */
    pm_start(pm);

    memcpy(pm->canvas, deltak_0, sizeof(pm->canvas[0]) * pm->allocsize);

    /* FIXME: Do it */
    /* paint to mesh */
    fastpm_particle_to_mesh(pm, pdata);

    pm_r2c(pm);

    /* copy out the results */
    memcpy(deltak_1, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);

    pm_stop(pm);
}

static int 
fastpm_particle_to_mesh(PM * pm, PMStore * p) 
{
    /* After this function, pm->canvas contains the real space density field */
    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
    };

    pm_append_ghosts(&pgd);

    pm_paint(pm, p, p->np + pgd.nghosts);

    pm_destroy_ghosts(&pgd);
}

