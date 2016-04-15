#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
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
fastpm_hmc_za_init(FastPMHMCZA * self, MPI_Comm comm)
{
    if(self->LPTOrder == 1)
        self->solver.USE_DX1_ONLY = 1;
    else
        self->solver.USE_DX1_ONLY = 0;

    double alloc_factor = 2.0;
    fastpm_2lpt_init(&self->solver, self->Nmesh, self->Ngrid, self->BoxSize, alloc_factor, comm);
    self->pm = self->solver.pm;
    self->delta_ic_k = pm_alloc(self->solver.pm);
    self->rho_final_x = pm_alloc(self->solver.pm);
}

void
fastpm_hmc_za_destroy(FastPMHMCZA * self)
{
    pm_free(self->solver.pm, self->rho_final_x);
    pm_free(self->solver.pm, self->delta_ic_k);
    fastpm_2lpt_destroy(&self->solver);
}

void 
fastpm_hmc_za_evolve(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic /* IC in k-space*/
    )
{
    FastPM2LPTSolver * solver = &self->solver;
    /* Evolve with ZA for HMC, smoothed by sml and deconvolve CIC */
    fastpm_2lpt_evolve(solver, delta_ic, 1.0, self->OmegaM);

    pm_assign(solver->pm, delta_ic, self->delta_ic_k);

    if(self->IncludeRSD) {
        fastpm_info("Using RSD along z\n");
        Cosmology c = {
                    .OmegaM = self->OmegaM,
                    .OmegaLambda = 1 - self->OmegaM,
        };

        double f1 = DLogGrowthFactor(1.0, c);
        double D1 = GrowthFactor(1.0, c);
        ptrdiff_t i;
        for(i = 0; i < solver->p->np; i ++) {
            solver->p->x[i][2] += f1 * D1 * solver->p->dx1[i][2];
        }
    }

    FastPMFloat * delta_final = self->rho_final_x;

    fastpm_utils_paint(solver->pm, solver->p, NULL, delta_final ,NULL, 0);

    if(self->SmoothingLength > 0)
        fastpm_apply_smoothing_transfer(solver->pm, delta_final, delta_final, self->SmoothingLength);
    if(self->KThreshold > 0)
        fastpm_apply_lowpass_transfer(solver->pm, delta_final, delta_final, self->KThreshold);
    if(self->DeconvolveCIC)
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
    FastPMFloat * sigma_x /* sigma_x in x-space*/
    )
{
    FastPM2LPTSolver * solver = &self->solver;
    FastPMFloat * model_x = self->rho_final_x;

    double chi2 = 0;

    PMXIter iter;
    for(pm_xiter_init(solver->pm, &iter);
       !pm_xiter_stop(&iter);
        pm_xiter_next(&iter)) {
        double diff = (model_x[iter.ind] - data_x[iter.ind]);
        diff *= diff;
        if(sigma_x)
            diff /= (sigma_x[iter.ind] * sigma_x[iter.ind]);
        chi2 += diff;
    }
    MPI_Allreduce(MPI_IN_PLACE, &chi2, 1, MPI_DOUBLE, MPI_SUM, self->solver.comm);
    chi2 /= pm_norm(solver->pm);
    return chi2;
}

void
fastpm_hmc_za_force(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * Fk    /* (out) hmc force in fourier space */
    )
{
    FastPM2LPTSolver * solver = &self->solver;

    FastPMFloat * model_x = self->rho_final_x;
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

    if(self->SmoothingLength > 0)
        fastpm_apply_smoothing_transfer(solver->pm, Fk, Fk, self->SmoothingLength);

    if(self->KThreshold > 0)
        fastpm_apply_lowpass_transfer(solver->pm, Fk, Fk, self->KThreshold);

    if(self->DeconvolveCIC)
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
            fastpm_info("Using RSD along z\n");
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
    }

    /* now we paint \Psi by the lagrangian position q */

    memset(Fk, 0, sizeof(Fk[0]) * pm_size(solver->pm));

    double fac = pm_norm(solver->pm) / pow(pm_boxsize(solver->pm)[0], 3);
    ptrdiff_t i;
    for(i = 0; i < solver->p->np; i ++) {
        for(d = 0; d < 3; d ++) {
        //    solver->p->q[i][d] += 0.25 * solver->boxsize / solver->nc;
        }
    }
    fastpm_utils_paint(solver->pm, solver->p, workspace2, NULL, get_lagrangian_position, 0);
    fastpm_utils_dump(solver->pm, "dump", workspace2);
    for(d = 0; d < 3; d ++) {
        fastpm_utils_paint(solver->pm, solver->p, NULL, workspace2, get_lagrangian_position, ACC[d]);
        fastpm_apply_za_hmc_force_transfer(solver->pm, workspace2, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm_size(solver->pm); ind ++) {
            /* Wang's magic factor of 2 in 1301.1348 is doubled because our chisq per ddof is approaching 1, not half.
             * We do not put it in in hmc_force_2lpt_transfer */
            Fk[ind] += 2 * 2 * fac * workspace[ind];
        }
    }
    for(i = 0; i < solver->p->np; i ++) {
        for(d = 0; d < 3; d ++) {
         //   solver->p->q[i][d] -= 0.25 * solver->boxsize / solver->nc;
        }
    }
    pm_free(solver->pm, workspace2);
    pm_free(solver->pm, workspace);
}

