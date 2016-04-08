#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/transfer.h>
#include <fastpm/hmc.h>

static void 
get_lagrangian_position(PMStore * p, ptrdiff_t index, double pos[3]) 
{
    pos[0] = p->q[index][0];
    pos[1] = p->q[index][1];
    pos[2] = p->q[index][2];
}
static void 
get_position(PMStore * p, ptrdiff_t index, double pos[3]) 
{
    pos[0] = p->x[index][0];
    pos[1] = p->x[index][1];
    pos[2] = p->x[index][2];
}

void 
fastpm_hmc_za_init(FastPMHMCZA * self, 
    int nc,
    double boxsize,
    double OmegaM, 
    MPI_Comm comm) 
{
    double alloc_factor = 2.0;
    self->OmegaM = OmegaM;
    self->sml = 0;
    self->kth = 0;
    self->solver.USE_DX1_ONLY = 1;
    fastpm_2lpt_init(&self->solver, nc, boxsize, alloc_factor, comm);
    self->decic = 0;
    self->pm = self->solver.pm;
}

void 
fastpm_hmc_za_destroy(FastPMHMCZA * self) 
{
    fastpm_2lpt_destroy(&self->solver);
}

void 
fastpm_hmc_za_evolve(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic, /* IC in k-space*/
    FastPMFloat * delta_final /* final in x-space*/
    ) 
{
    FastPM2LPTSolver * solver = &self->solver;
    /* Evolve with ZA for HMC, smoothed by sml and deconvolve CIC */
    fastpm_2lpt_evolve(solver, delta_ic, 1.0, self->OmegaM);
    fastpm_utils_paint(solver->pm, solver->p, NULL, delta_final, NULL, 0);

    if(self->sml > 0)
        fastpm_apply_smoothing_transfer(solver->pm, delta_final, delta_final, self->sml);
    if(self->kth > 0)
        fastpm_apply_lowpass_transfer(solver->pm, delta_final, delta_final, self->kth);
    if(self->decic)
        fastpm_apply_decic_transfer(solver->pm, delta_final, delta_final);

    pm_c2r(solver->pm, delta_final);
    ptrdiff_t ind;
    /* inv volume of a cell, to convert to density */
    double fac = (pm_norm(solver->pm) / pow(pm_boxsize(solver->pm)[0], 3));
    for(ind = 0; ind < pm_size(solver->pm); ind++) {
        delta_final[ind] *= fac;
    }
}

double
fastpm_hmc_za_chisq(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * model_x, /* rhop in x-space*/
    FastPMFloat * sigma_x /* sigma_x in x-space*/
    )
{
    FastPM2LPTSolver * solver = &self->solver;
    double chi2 = 0;

    PMXIter iter;
    for(pm_xiter_init(solver->pm, &iter);
       !pm_xiter_stop(&iter);
        pm_xiter_next(&iter)) {
        double diff = (model_x[iter.ind] - data_x[iter.ind]);
        diff *= diff;
        if(sigma_x)
            diff /= (2 * sigma_x[iter.ind] * sigma_x[iter.ind]);
        chi2 += diff;
        //printf("%td\n", iter.ind);
    }
    MPI_Allreduce(MPI_IN_PLACE, &chi2, 1, MPI_DOUBLE, MPI_SUM, self->solver.comm);
    chi2 /= pm_norm(solver->pm);
    return chi2;
}

void
fastpm_hmc_za_force(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * model_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * Fk    /* (out) hmc force in fourier space */
    )
{
    FastPM2LPTSolver * solver = &self->solver;
    int d;

    FastPMFloat * workspace = pm_alloc(solver->pm);
    FastPMFloat * workspace2 = pm_alloc(solver->pm);

    ptrdiff_t ind;
    for(ind = 0; ind < pm_size(solver->pm); ind ++) {
        workspace[ind] = model_x[ind] - data_x[ind];
        if(sigma_x)
            workspace[ind] /= sigma_x[ind] * sigma_x[ind];
    }

    pm_r2c(solver->pm, workspace, Fk);

    if(self->sml > 0)
        fastpm_apply_smoothing_transfer(solver->pm, Fk, Fk, self->sml);
    if(self->kth > 0)
        fastpm_apply_lowpass_transfer(solver->pm, Fk, Fk, self->kth);

    if(self->decic)
        fastpm_apply_decic_transfer(solver->pm, Fk, Fk);

    /* Fk contains rhod_k at this point */

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        fastpm_apply_diff_transfer(solver->pm, Fk, workspace, d);

        /* workspace stores \Gamma(k) = i k \rho_d */

        pm_c2r(solver->pm, workspace);

        fastpm_utils_readout(solver->pm, solver->p, workspace, get_position, ACC[d]);

        /*FIXME: Add RSD f */
        if(self->IncludeRSD && d == 2) {
            Cosmology c = {
                        .OmegaM = self->OmegaM,
                        .OmegaLambda = 1 - self->OmegaM,
            };

            double f1 = DLogGrowthFactor(1.0, c);
            double D1 = GrowthFactor(1.0, c);
            ptrdiff_t i;
            for(i = 0; i < solver->p->np; i ++) {
                solver->p->acc[i][d] *= (1 + f1 * D1);
            }
        }
        fastpm_info("d= %d ACC = %d\n", d, ACC[d]);
    }

    /* now we paint \Psi by the lagrangian position q */

    memset(Fk, 0, sizeof(Fk[0]) * pm_size(solver->pm));

    double fac = pm_norm(solver->pm) / pow(pm_boxsize(solver->pm)[0], 3);

    for(d = 0; d < 3; d ++) {
        fastpm_utils_paint(solver->pm, solver->p, NULL, workspace2, get_lagrangian_position, ACC[d]);
        fastpm_apply_za_hmc_force_transfer(solver->pm, workspace2, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm_size(solver->pm); ind ++) {
            /* Wang's magic factor of 2 in 1301.1348.
             * We do not put it in in hmc_force_2lpt_transfer */
            Fk[ind] += 2 * fac * workspace[ind];
        }
    }
    pm_free(solver->pm, workspace2);
    pm_free(solver->pm, workspace);
}

