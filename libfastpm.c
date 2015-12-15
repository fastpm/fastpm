/* libfastpm: */
#include <mpi.h>
#include <string.h>
#include "libfastpm.h"
#include "msg.h"
#include "parameters.h"
#include "pmsteps.h"
#include "pmic.h"

static int 
fastpm_particle_to_mesh(PM * pm, PMStore * p);

int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm) 
{
    msg_init(comm);
    msg_set_loglevel(verbose);

    pm_store_init(p);

    pm_store_alloc_evenly(p, pow(nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_ACC | PACK_DX1 | PACK_DX2 | PACK_Q,
        2.0, comm);

    return 0;
}

void
fastpm_init_pm(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm) 
{
    pm_init_simple(pm, p, Ngrid, BoxSize, comm);
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
fastpm_fill_deltak(PM * pm, float_t * deltak, int seed, fastpm_pkfunc pk, void * pkdata) 
{
    pm_ic_fill_gaussian_gadget(pm, seed, pk, pkdata);
    memcpy(deltak, pm->canvas, sizeof(pm->canvas[0]) * pm->allocsize);
}

void 
fastpm_evolve_2lpt(PM * pm, PMStore * pdata, 
        double a, double omega_m, 
        float_t * deltak_0)
{
    /* evolve particles by 2lpt to time a. pm->canvas contains rho(x, a) */
    double shift[3] = {0, 0, 0};

    pm_store_set_lagrangian_position(pdata, pm, shift);

    memcpy(pm->canvas, deltak_0, sizeof(pm->canvas[0]) * pm->allocsize);

    pm_2lpt_main(pm, pdata, shift);

    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */
    int i;

    /* predict particle positions by 2lpt */
    pm_2lpt_evolve(a, pdata, omega_m);

    fastpm_particle_to_mesh(pm, pdata);
}

static int 
to_rank(void * pdata, ptrdiff_t i, void * data) 
{
    PMStore * p = (PMStore *) pdata;
    PM * pm = (PM*) data;
    double pos[3];
    p->iface.get_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
}

void 
fastpm_evolve_pm(PM * basepm, VPM * vpm_list, 
        PMStore * pdata, 
        double time_step[],
        int n_time_step,
        double omega_m, 
        float_t * deltak_0, float_t * deltak_1, MPI_Comm comm) 
{
    double shift[3] = {0, 0, 0};

    pm_store_set_lagrangian_position(pdata, basepm, shift);

    pm_start(basepm);

    memcpy(basepm->canvas, deltak_0, sizeof(basepm->canvas[0]) * basepm->allocsize);

    pm_2lpt_main(basepm, pdata, shift);

    /* predict particle positions by 2lpt */
    pm_2lpt_evolve(time_step[0], pdata, omega_m);

    PMStepper stepper;

    stepping_init(&stepper, omega_m, FORCE_MODE_PM, 1);

    int istep;
    int nsteps = n_time_step;

    for (istep = 0; istep < nsteps; istep++) {
        double a_v, a_x, a_v1, a_x1;

        pm_get_times(istep, time_step, n_time_step,
            &a_x, &a_x1, &a_v, &a_v1);

        /* Find the Particle Mesh to use for this time step */
        VPM * vpm = vpm_find(vpm_list, a_x);
        PM * pm = &vpm->pm;

        /* watch out: boost the density since mesh is finer than grid */
        double density_factor =  pow(vpm->pm_nc_factor, 3); 

        pm_store_wrap(pdata, pm->BoxSize);

        pm_store_decompose(pdata, to_rank, pm, comm);

        pm_start(pm);

        pm_calculate_forces(pdata, pm, density_factor);

        pm_stop(pm);

        stepping_kick(&stepper, pdata, pdata, a_v, a_v1, a_x);

        stepping_drift(&stepper, pdata, pdata, a_x, a_x1, a_v1);
    }

    /* paint to mesh */
    fastpm_particle_to_mesh(basepm, pdata);

    pm_r2c(basepm);

    /* copy out the results */
    memcpy(deltak_1, basepm->canvas, sizeof(basepm->canvas[0]) * basepm->allocsize);

    pm_stop(basepm);
}


static void 
get_lagrangian_position(void * pdata, ptrdiff_t index, double pos[3]) 
{
    PMStore * p = (PMStore *)pdata;
    pos[0] = p->q[index][0];
    pos[1] = p->q[index][1];
    pos[2] = p->q[index][2];
}

void fastpm_apply_diff_transfer(PM * pm, int dir) {

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k;

            /* - i k[d] */
            pm->workspace[ind + 0] =   pm->canvas[ind + 1] * (k_finite);
            pm->workspace[ind + 1] = - pm->canvas[ind + 0] * (k_finite);
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);

}

void fastpm_apply_hmc_force_2lpt_transfer(PM * pm, int dir) {

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k;
            double kk_finite = 0.;
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk;
            }

            if(kk_finite == 0)
            {
                pm->workspace[ind + 0] = 0;
                pm->workspace[ind + 1] = 0;
            }
            else
            {
                /* - i k[d] / k**2 */
                pm->workspace[ind + 0] =   pm->canvas[ind + 1] * (k_finite / kk_finite);
                pm->workspace[ind + 1] = - pm->canvas[ind + 0] * (k_finite / kk_finite);
            }
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);

}

void 
fastpm_derivative_2lpt(PM * pm, 
        PMStore * p, /* Current position (x) saved in -> x */
        float_t * rhop_x, /* rhop in x-space*/
        float_t * Fk     /* (out) hmc force in fourier space */
        )
{
    int d;

    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
        .get_position = p->iface.get_position,
        .nghosts = 0,
    };

    pm_append_ghosts(&pgd);

    pm_paint(pm, p, p->np + pgd.nghosts);

    ptrdiff_t ind;

    for(ind = 0; ind < pm->allocsize; ind ++) {
        pm->workspace[ind] -= rhop_x[ind];
    }

    pm_r2c(pm);

    /* canvas contains rhod_k at this point */

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};

    for(d = 0; d < 3; d ++) {
        fastpm_apply_diff_transfer(pm, d);

        /* workspace stores \Gamma(k) = -i k \rho_d */

        pm_c2r(pm);
        
        int i;
        /* acc stores \Gamma(x) := \Psi(q) */
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            p->acc[i][d] = pm_readout_one(pm, p, i) / pm->Norm;
        }
        pm_reduce_ghosts(&pgd, ACC[d]);
    }

    pm_destroy_ghosts(&pgd);

    /* now we paint \Psi by the lagrangian position q */
    pgd.get_position = get_lagrangian_position;
    pgd.attributes = PACK_Q;

    pm_append_ghosts(&pgd);

    memset(pm->workspace, 0, sizeof(pm->workspace[0]) * pm->allocsize);
    memset(Fk, 0, sizeof(Fk[0]) * pm->allocsize);

    for(d = 0; d < 3; d ++) {
        int i;
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            double pos[3];
            get_lagrangian_position(p, i, pos);
            pm_paint_pos(pm, pos, p->acc[i][d]);
        }
        pm_r2c(pm);
        fastpm_apply_hmc_force_2lpt_transfer(pm, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm->allocsize; ind ++) {
            Fk[ind] += 2 * pm->workspace[ind] / pm->Norm; 
            /*Wang's magic factor of 2 in 1301.1348. 
             * We do not put it in in hmc_force_2lpt_transfer */
        }
    }
    pm_destroy_ghosts(&pgd);
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
        .get_position = p->iface.get_position,
        .nghosts = 0,
    };

    pm_append_ghosts(&pgd);

    pm_paint(pm, p, p->np + pgd.nghosts);

    pm_destroy_ghosts(&pgd);
}

