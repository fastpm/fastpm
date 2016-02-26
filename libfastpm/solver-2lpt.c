/* libfastpm: */
#include <mpi.h>
#include <string.h>

#include <fastpm/libfastpm.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "transfer.h"

int 
fastpm_2lpt_init(FastPM2LPTSolver * solver, int nc, double boxsize, double alloc_factor, MPI_Comm comm)
{
    solver->p = malloc(sizeof(PMStore));
    solver->pm = malloc(sizeof(PM));
    pm_store_init(solver->p);

    pm_store_alloc_evenly(solver->p, pow(nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_ACC | PACK_DX1 | PACK_DX2 | PACK_Q,
        alloc_factor, comm);

    pm_init_simple(solver->pm, solver->p, nc, boxsize, comm);
    solver->boxsize = boxsize;
    solver->nc = nc;
    solver->comm = comm;
    return 0;
}

void
fastpm_2lpt_destroy(FastPM2LPTSolver * solver) 
{
    pm_destroy(solver->pm);
    pm_store_destroy(solver->p);
    free(solver->pm);
    free(solver->p);   
}

void 
fastpm_2lpt_evolve(FastPM2LPTSolver * solver, 
        FastPMFloat * delta_k_i, double aout, double omega_m)
{
    /* evolve particles by 2lpt to time a. pm->canvas contains rho(x, a) */
    double shift[3] = {0, 0, 0
            //solver->boxsize / solver->nc * 0.5,
            //solver->boxsize / solver->nc * 0.5,
            //solver->boxsize / solver->nc * 0.5
    };

    pm_store_set_lagrangian_position(solver->p, solver->pm, shift);

    pm_2lpt_solve(solver->pm, delta_k_i, solver->p, shift);

    pm_store_summary(solver->p, solver->pm->Comm2D);

    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */

    /* predict particle positions by 2lpt */
    pm_2lpt_evolve(aout, solver->p, omega_m, solver->USE_DX1_ONLY);
}

static void 
get_lagrangian_position(void * pdata, ptrdiff_t index, double pos[3]) 
{
    PMStore * p = (PMStore *)pdata;
    pos[0] = p->q[index][0];
    pos[1] = p->q[index][1];
    pos[2] = p->q[index][2];
}

void fastpm_apply_hmc_force_2lpt_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) {

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double k_finite = kiter.fac[dir][kiter.iabs[dir]].k_finite;
            double kk_finite = 0.;
            double cic = 1.0;            
            for(d = 0; d < 3; d++) {
                kk_finite += kiter.fac[d][kiter.iabs[d]].kk_finite;
                /*  cic *= kiter.fac[d][kiter.iabs[d]].cic; */
            }
            if(kk_finite == 0)
            {
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            }
            else
            {
                /* - i k[d] / k**2 */
                to[kiter.ind + 0] =   from[kiter.ind + 1] * (k_finite / kk_finite / cic);
                to[kiter.ind + 1] = - from[kiter.ind + 0] * (k_finite / kk_finite / cic);
            }
        }
    }
}

void 
fastpm_2lpt_hmc_force(FastPM2LPTSolver * solver,
        FastPMFloat * data_x, /* rhop in x-space*/
        FastPMFloat * sigma_x, /* sigma_x in x-space*/
        FastPMFloat * Fk,    /* (out) hmc force in fourier space */
        double sml
        )
{
    int d;

    FastPMFloat * workspace = pm_alloc(solver->pm);
    FastPMFloat * workspace2 = pm_alloc(solver->pm);

    PMGhostData * pgd = pm_ghosts_create(solver->pm, solver->p, PACK_POS, NULL);

    /* Note the 1.0 see the other comment on pm_paint in this file.*/
    pm_paint(solver->pm, workspace, solver->p, solver->p->np + pgd->nghosts, 1.0);

    ptrdiff_t ind;

    pm_r2c(solver->pm, workspace, Fk);

    fastpm_apply_smoothing_transfer(solver->pm, Fk, workspace, sml);

    pm_c2r(solver->pm, workspace);

    for(ind = 0; ind < solver->pm->allocsize; ind ++) {
        workspace[ind] /= pm_norm(solver->pm);
        workspace[ind] -= data_x[ind];
        if(sigma_x)
            workspace[ind] /= sigma_x[ind] * sigma_x[ind];
    }

    pm_r2c(solver->pm, workspace, Fk);

    /* Fk contains rhod_k at this point */

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};

    for(d = 0; d < 3; d ++) {
        fastpm_apply_diff_transfer(solver->pm, Fk, workspace, d);

        /* workspace stores \Gamma(k) = -i k \rho_d */

        pm_c2r(solver->pm, workspace);
        
        int i;
        /* acc stores \Gamma(x) := \Psi(q) */
        for(i = 0; i < solver->p->np + pgd->nghosts; i ++) {
            solver->p->acc[i][d] = pm_readout_one(solver->pm, workspace, solver->p, i) / solver->pm->Norm;
        }
        pm_ghosts_reduce(pgd, ACC[d]);
    }

    pm_ghosts_free(pgd);

    /* now we paint \Psi by the lagrangian position q */

    pgd = pm_ghosts_create(solver->pm, solver->p, PACK_Q | PACK_ACC, get_lagrangian_position);

    memset(Fk, 0, sizeof(Fk[0]) * solver->pm->allocsize);

    for(d = 0; d < 3; d ++) {
        memset(workspace, 0, sizeof(workspace[0]) * solver->pm->allocsize);
        memset(workspace2, 0, sizeof(workspace2[0]) * solver->pm->allocsize);
        int i;
        for(i = 0; i < solver->p->np + pgd->nghosts; i ++) {
            double pos[3];
            get_lagrangian_position(solver->p, i, pos);
            pm_paint_pos(solver->pm, workspace, pos, solver->p->acc[i][d]);
        }
        pm_r2c(solver->pm, workspace, workspace2);
        fastpm_apply_hmc_force_2lpt_transfer(solver->pm, workspace2, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < solver->pm->allocsize; ind ++) {
            Fk[ind] += - 2 * workspace[ind] / solver->pm->Norm; 
            /* Wang's magic factor of 2 in 1301.1348. 
             * We do not put it in in hmc_force_2lpt_transfer */
        }
    }
    pm_ghosts_free(pgd);
    pm_free(solver->pm, workspace2);
    pm_free(solver->pm, workspace);
}

