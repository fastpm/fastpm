#include <string.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

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

static void
fastpm_pm_init(FastPM * fastpm, PMStore * p, int nmesh, int nc, double boxsize, double alloc_factor, double omega_m, MPI_Comm comm)
{
    fastpm->nc = nc;
    fastpm->boxsize = boxsize;
    fastpm->alloc_factor = alloc_factor;
    fastpm->omega_m = omega_m;
    fastpm->vpminit = (VPMInit[]) {
        {.a_start = 0, .pm_nc_factor = 1},
        {.a_start = -1, .pm_nc_factor = 0},
    };
    fastpm->FORCE_TYPE = FASTPM_FORCE_FASTPM;
    fastpm->USE_NONSTDDA = 0;
    fastpm->USE_EXTERNAL_PSTORE = p;
    fastpm->USE_MODEL = FASTPM_MODEL_NONE;
    fastpm->USE_SHIFT = 0;
    fastpm->nLPT = 2.5;
    fastpm->K_LINEAR = 4;
    fastpm->KERNEL_TYPE = FASTPM_KERNEL_EASTWOOD;
    fastpm_init(fastpm, 0, 0, comm);
}

void 
fastpm_hmc_za_init(FastPMHMCZA * self, MPI_Comm comm)
{
    switch(self->LPTOrder) {
        case 1:
            self->solver.USE_DX1_ONLY = 1;
            self->pm_solver.USE_DX1_ONLY = 1;
        break;
        case 2:
            self->solver.USE_DX1_ONLY = 0;
            self->pm_solver.USE_DX1_ONLY = 0;
        break;
        default:
            fastpm_raise(-1, "Wrong LPT Order, can only be 1 or 2\n");
    }

    fastpm_2lpt_init(&self->solver, self->Nmesh, self->Ngrid, self->BoxSize, self->AllocFactor, comm);
    self->solver.USE_SHIFT = 0;

    fastpm_pm_init(&self->pm_solver, self->solver.p, self->Nmesh, self->Ngrid, self->BoxSize, self->AllocFactor, self->OmegaM, comm);

    /* FIXME: create a new pm object */
    self->pm = self->solver.pm;

    /* We will set p after evolve is called. p contains the displacement field! */
    self->p = NULL;

    self->delta_ic_k = pm_alloc(self->solver.pm);
    self->rho_final_x = pm_alloc(self->solver.pm);
    self->transfer_function = pm_alloc(self->solver.pm);
    if(self->IncludeRSD) {
        fastpm_info("Using RSD along z\n");
    }
}


void
fastpm_hmc_za_destroy(FastPMHMCZA * self)
{
    pm_free(self->solver.pm, self->transfer_function);
    pm_free(self->solver.pm, self->rho_final_x);
    pm_free(self->solver.pm, self->delta_ic_k);
    fastpm_destroy(&self->pm_solver);
    fastpm_2lpt_destroy(&self->solver);
}

static void
fastpm_hmc_za_evolve_internal(
    FastPMHMCZA * self,
    int Nsteps, /* 0 for LPT, > 0 for PM */
    FastPMFloat * delta_final
    )
{
    /* This is the pure N-body dynamics, no smoothing etc */
    /* First run a PT simulation*/
    if(Nsteps <= 0) {
        FastPM2LPTSolver * solver = &self->solver;
        fastpm_2lpt_evolve(solver, self->delta_ic_k, 1.0, self->OmegaM);

        self->p = solver->p;
    } else {
        FastPM * solver = &self->pm_solver;
        double time_step[Nsteps];
        double ainit = 0.1;
        double afinal = 1.0;
        int i;
        for(i = 0; i < Nsteps; i ++) {
            time_step[i] = ainit + 1.0 * i / (Nsteps - 1) * (afinal - ainit);
        }
        time_step[Nsteps - 1] = afinal;
        fastpm_setup_ic(solver, self->delta_ic_k);
        fastpm_evolve(solver, time_step, Nsteps);

        self->p = solver->p;
    }
    if(self->IncludeRSD) {
        Cosmology c = {
            .OmegaM = self->OmegaM,
            .OmegaLambda = 1 - self->OmegaM,
        };
        double RSD = 1.0 / (self->p->a_x * self->p->a_x * HubbleEa(self->p->a_x, c));
        ptrdiff_t i;
        for(i = 0; i < self->p->np; i ++) {
            self->p->x[i][2] += self->p->v[i][2] * RSD;
        }
    }

    fastpm_utils_paint(self->pm, self->p, NULL, delta_final, NULL, 0);
}

void
fastpm_hmc_za_calibrate_transferfunction(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic, /* IC in k-space*/
    int Nsteps,
    FastPMPowerSpectrum * transfer)
{
    FastPMFloat * deltalpt = pm_alloc(self->pm);
    FastPMFloat * deltapm = pm_alloc(self->pm);

    pm_assign(self->pm, delta_ic, self->delta_ic_k);

    fastpm_hmc_za_evolve_internal(self, Nsteps, deltalpt);

    pm_assign(self->pm, delta_ic, self->delta_ic_k);

    fastpm_hmc_za_evolve_internal(self, 0, deltapm);
    fastpm_transferfunction_init(transfer, self->pm, deltalpt, deltapm);
}

void
fastpm_hmc_za_evolve(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic, /* IC in k-space*/
    int Nsteps /* 0 for LPT, > 0 for PM */
    )
{
    pm_assign(self->pm, delta_ic, self->delta_ic_k);

    FastPMFloat * delta_final = self->rho_final_x;

    fastpm_hmc_za_evolve_internal(self, Nsteps, delta_final);

    if(self->SmoothingLength > 0)
        fastpm_apply_smoothing_transfer(self->pm, delta_final, delta_final, self->SmoothingLength);
    if(self->KThreshold > 0)
        fastpm_apply_lowpass_transfer(self->pm, delta_final, delta_final, self->KThreshold);
    if(self->DeconvolveCIC)
        fastpm_apply_decic_transfer(self->pm, delta_final, delta_final);
    if(self->TransferFunction.func)
        fastpm_apply_any_transfer(self->pm, delta_final, delta_final, self->TransferFunction.func, self->TransferFunction.data);

    pm_c2r(self->pm, delta_final);

    //  inv volume of a cell, to convert to density
    double fac = (pm_norm(self->pm) / pow(pm_boxsize(self->pm)[0], 3));
    ptrdiff_t ind;
    for(ind = 0; ind < pm_size(self->pm); ind++) {
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
    FastPMFloat * model_x = self->rho_final_x;

    double chi2 = 0;

    PMXIter iter;
    for(pm_xiter_init(self->pm, &iter);
       !pm_xiter_stop(&iter);
        pm_xiter_next(&iter)) {
        double diff = (model_x[iter.ind] - data_x[iter.ind]);
        diff *= diff;
        if(sigma_x)
            diff /= (sigma_x[iter.ind] * sigma_x[iter.ind]);
        chi2 += diff;
    }
    MPI_Allreduce(MPI_IN_PLACE, &chi2, 1, MPI_DOUBLE, MPI_SUM, self->solver.comm);
    chi2 /= pm_norm(self->pm);
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
    ptrdiff_t ind;

    FastPMFloat * rhodk = pm_alloc(self->pm);
    FastPMFloat * workspace = pm_alloc(self->pm);

    fastpm_hmc_za_force_rhodk(self, data_x, sigma_x, rhodk);

    memset(Fk, 0, sizeof(Fk[0]) * pm_size(self->pm));

    /* First order */
    if(self->LPTOrder >= 1) {
        fastpm_hmc_za_force_s1(self, rhodk, workspace);

        for(ind = 0; ind < pm_size(self->pm); ind++) {
            Fk[ind] += workspace[ind];
        }
    }
    /* Second order */
    if(self->LPTOrder >= 2) {
        fastpm_hmc_za_force_s2(self, Fk, workspace);

        for(ind = 0; ind < pm_size(self->pm); ind++) {
            Fk[ind] += workspace[ind];
        }
    }

    pm_free(self->pm, workspace);
    pm_free(self->pm, rhodk);
}

void
fastpm_hmc_za_force_rhodk(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * rhodk    /* (out) rhodk in fourier space */
    )
{

    FastPMFloat * model_x = self->rho_final_x;

    FastPMFloat * workspace = pm_alloc(self->pm);

    ptrdiff_t ind;
    for(ind = 0; ind < pm_size(self->pm); ind ++) {
        workspace[ind] = model_x[ind] - data_x[ind];
        if(sigma_x)
            workspace[ind] /= sigma_x[ind] * sigma_x[ind];
    }

    pm_r2c(self->pm, workspace, rhodk);

    pm_free(self->pm, workspace);

    if(self->SmoothingLength > 0)
        fastpm_apply_smoothing_transfer(self->pm, rhodk, rhodk, self->SmoothingLength);

    if(self->KThreshold > 0)
        fastpm_apply_lowpass_transfer(self->pm, rhodk, rhodk, self->KThreshold);

    if(self->DeconvolveCIC)
        fastpm_apply_decic_transfer(self->pm, rhodk, rhodk);

    if(self->TransferFunction.func)
        fastpm_apply_any_transfer(self->pm, rhodk, rhodk, self->TransferFunction.func, self->TransferFunction.data);

}
void
fastpm_hmc_za_force_s1(
    FastPMHMCZA * self,
    FastPMFloat * rhodk,
    FastPMFloat * Fk1    /* (out) hmc force for s2 in fourier space */
    )
{
    int d;

    FastPMFloat * workspace = pm_alloc(self->pm);

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        fastpm_apply_diff_transfer(self->pm, rhodk, workspace, d);

        /* workspace stores \Gamma(k) = i k \rho_d */

        pm_c2r(self->pm, workspace);

        fastpm_utils_readout(self->pm, self->p, workspace, get_position, ACC[d]);

        /*FIXME: Add RSD f */
        if(self->IncludeRSD && d == 2) {
            Cosmology c = {
                        .OmegaM = self->OmegaM,
                        .OmegaLambda = 1 - self->OmegaM,
            };

            double f1 = DLogGrowthFactor(1.0, c);
            double D1 = GrowthFactor(1.0, c);
            ptrdiff_t i;
            for(i = 0; i < self->p->np; i ++) {
                self->p->acc[i][d] *= (1 + f1 * D1);
            }
        }
    }

    /* now we paint \Psi by the lagrangian position q */

    memset(Fk1, 0, sizeof(Fk1[0]) * pm_size(self->pm));

    double fac = pm_norm(self->pm) / pow(pm_boxsize(self->pm)[0], 3);

    for(d = 0; d < 3; d ++) {
        fastpm_utils_paint(self->pm, self->p, NULL, workspace, get_lagrangian_position, ACC[d]);
        fastpm_apply_laplace_transfer(self->pm, workspace, workspace);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            /* Wang's magic factor of 2 in 1301.1348 is doubled because our chisq per ddof is approaching 1, not half.
             * We do not put it in in hmc_force_2lpt_transfer */

            /* negative sign because the force is - grad */
            Fk1[ind] += - 2 * 2 * fac * workspace[ind];
        }
    }
    pm_free(self->pm, workspace);
}

void
fastpm_hmc_za_force_s2(
    FastPMHMCZA * self,
    FastPMFloat * Fk1,   /* (in) hmc force for s1 in fourier space */
    FastPMFloat * Fk2    /* (out) hmc force for s2 in fourier space */
    )
{
    /* This function is originally written by Chirag Modi. <chirag@berkeley.edu> */

    int d;
    ptrdiff_t ind;


    FastPMFloat * Fpsi = pm_alloc(self->pm);
    FastPMFloat * source = pm_alloc(self->pm);
    FastPMFloat * workspace = pm_alloc(self->pm);

    pm_assign(self->pm, Fk1, Fpsi);
    for(ind = 0; ind < pm_size(self->pm); ind ++) {
        /* FIXME: put in D2 */
        Fpsi[ind] *= 3.0 / 7.0;
    }

    /* 2LPT derivative begins here , carrying over Fk2 from above */

    pm_c2r(self->pm, Fpsi);

    memset(Fk2, 0, sizeof(Fk2[0]) * pm_size(self->pm));

    /* diagonal elements */
    for(d = 0; d < 3; d++){

        fastpm_apply_laplace_transfer(self->pm, self->delta_ic_k, workspace);

        fastpm_apply_diff_transfer(self->pm, workspace, workspace, d);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, d);
        pm_c2r(self->pm, workspace);

        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            workspace[ind]  = Fpsi[ind] * workspace[ind];
        }

        pm_r2c(self->pm, workspace, source);

        fastpm_apply_laplace_transfer(self->pm, source, source);


        fastpm_apply_diff_transfer(self->pm, source, workspace, (d+1)%3);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, (d+1)%3);

        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            Fk2[ind] += workspace[ind];
        }

        fastpm_apply_diff_transfer(self->pm, source, workspace, (d+2)%3);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, (d+2)%3);

        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            Fk2[ind] += workspace[ind];

        }
    }

    /* off - diagonal elements */
    for(d = 0; d < 3; d++){

        fastpm_apply_laplace_transfer(self->pm, self->delta_ic_k, workspace);

        fastpm_apply_diff_transfer(self->pm, workspace, workspace, (d+1)%3);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, (d+2)%3);

        pm_c2r(self->pm, workspace);

        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            workspace[ind]  = Fpsi[ind] * workspace[ind];
        }

        pm_r2c(self->pm, workspace, source);

        fastpm_apply_laplace_transfer(self->pm, source, source);


        fastpm_apply_diff_transfer(self->pm, source, workspace, (d+1)%3);
        fastpm_apply_diff_transfer(self->pm, workspace, workspace, (d+2)%3);

        for(ind = 0; ind < pm_size(self->pm); ind ++) {
            Fk2[ind] -= 2 * workspace[ind];
        }
    }

    pm_free(self->pm, workspace);
    pm_free(self->pm, source);
    pm_free(self->pm, Fpsi);
}














//////////////////////////   TRASH SAVED HERE ///////////////////////////////////


//From hmc_za_evolve, after declaring delta_final
    /*
    fastpm_utils_paint(solver->pm, solver->p, NULL, delta_final ,NULL, 0);

    if(self->SmoothingLength > 0)
        fastpm_apply_smoothing_transfer(solver->pm, delta_final, delta_final, self->SmoothingLength);
    if(self->KThreshold > 0)
        fastpm_apply_lowpass_transfer(solver->pm, delta_final, delta_final, self->KThreshold);
    if(self->DeconvolveCIC)
        fastpm_apply_decic_transfer(solver->pm, delta_final, delta_final);

    pm_c2r(solver->pm, delta_final);
    ptrdiff_t ind;
    //  inv volume of a cell, to convert to density
    double fac = (pm_norm(solver->pm) / pow(pm_boxsize(solver->pm)[0], 3));
    for(ind = 0; ind < pm_size(solver->pm); ind++) {
        delta_final[ind] *= fac;
    }
    */
    

    /*
    FastPM2LPTSolver * solver = &self->solver;
    // Evolve with ZA for HMC, smoothed by sml and deconvolve CIC 
    fastpm_2lpt_evolve(solver, delta_ic, 1.0, self->OmegaM);

    pm_assign(solver->pm, delta_ic, self->delta_ic_k);

    pmt = solver->pm;
    pt = solver->p;
    */
    
