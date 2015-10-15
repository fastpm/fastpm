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

    pm_store_alloc_evenly(p, pow(nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_ACC | PACK_DX1 | PACK_DX2 | PACK_Q,
        2.0, comm);

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

    int i;
    /* copy the lagrangian coordinates. */
    for(i = 0; i < pdata->np; i ++) {
        for(d = 0; d < 3; d ++) {
            pdata->q[i][d] = pdata->x[i][d] + shift[d];
        }
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
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k_finite;

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
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k_finite;
            double kk_finite = 0.;
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
            }

            /* - i k[d] / k**2 */
            pm->workspace[ind + 0] =   pm->canvas[ind + 1] * (k_finite / kk_finite);
            pm->workspace[ind + 1] = - pm->canvas[ind + 0] * (k_finite / kk_finite);
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);

}

void 
fastpm_derivative_2lpt(PM * pm, 
        PMStore * p, /* Current position (x) saved in -> x */
        real_t * rhod_k, /* rhod in fourier space */
        real_t * Fk,     /* (out) hmc force in fourier space */
        MPI_Comm comm) 
{
    int d;
    pm_start(pm);

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

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};

    for(d = 0; d < 3; d ++) {
        memcpy(pm->canvas, rhod_k, sizeof(pm->canvas[0]) * pm->allocsize);

        fastpm_apply_diff_transfer(pm, d);

        /* Canvas stores \Gamma(k) = -i k \rho_d */

        pm_c2r(pm);
        
        int i;
        /* acc stores \Gamma(x) := \Psi(q) */
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            p->acc[i][d] = pm_readout_one(pm, p, i);
        }
        pm_reduce_ghosts(&pgd, ACC[d]);
    }

    pm_destroy_ghosts(&pgd);

    /* now we paint \Psi by the lagrangian position q */
    pgd.get_position = get_lagrangian_position;
    pgd.attributes = PACK_Q;

    pm_append_ghosts(&pgd);

    memset(pm->canvas, 0, sizeof(pm->canvas[0]) * pm->allocsize);
    memset(Fk, 0, sizeof(Fk[0]) * pm->allocsize);

    for(d = 0; d < 3; d ++) {
        int i;
        for(i = 0; i < p->np + pgd.nghosts; p ++) {
            double pos[3];
            get_lagrangian_position(p, i, pos);
            pm_paint_pos(pm, pos, p->acc[i][d]);
        }
        pm_r2c(pm);
        fastpm_apply_hmc_force_2lpt_transfer(pm, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm->allocsize; ind ++) {
            Fk[ind] += pm->workspace[ind];
        }
    }
    pm_destroy_ghosts(&pgd);

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
        .get_position = p->iface.get_position,
        .nghosts = 0,
    };

    pm_append_ghosts(&pgd);

    pm_paint(pm, p, p->np + pgd.nghosts);

    pm_destroy_ghosts(&pgd);
}

