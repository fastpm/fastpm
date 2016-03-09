#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/transfer.h>
typedef struct {
    FastPM2LPTSolver solver;
    double OmegaM;
    double sml;
    int decic;
    PM * pm;
} FastPMHMCZA;

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
    double sml,
    MPI_Comm comm) 
{
    double alloc_factor = 2.0;
    self->OmegaM = OmegaM;
    self->sml = sml;
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

    FastPMFloat * tmp = pm_alloc(solver->pm);

    pm_r2c(solver->pm, delta_final, tmp);
    fastpm_apply_smoothing_transfer(solver->pm, tmp, delta_final, self->sml);

    pm_free(solver->pm, tmp);

    if(self->decic)
        fastpm_apply_decic_transfer(solver->pm, delta_final, delta_final);

    pm_c2r(solver->pm, delta_final);
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

    PMKIter kiter;
    for(pm_kiter_init(solver->pm, &kiter); 
        !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        ptrdiff_t ind = kiter.ind;
        double w = 1.0;
        if(kiter.iabs[0] == 0 ||
           kiter.iabs[1] == 0 ||
           kiter.iabs[2] == 0) {
            w = 0.5;
        }
        double diff = w * (model_x[ind] - data_x[ind]);
        diff *= diff;
        if(sigma_x)
            diff /= (2 * sigma_x[ind] * sigma_x[ind]);
        chi2 += diff;
    }

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

    fastpm_apply_smoothing_transfer(solver->pm, Fk, Fk, self->sml);

    if(self->decic)
        fastpm_apply_decic_transfer(solver->pm, Fk, Fk);

    /* Fk contains rhod_k at this point */

    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) { 
        fastpm_apply_diff_transfer(solver->pm, Fk, workspace, d);

        /* workspace stores \Gamma(k) = i k \rho_d */

        pm_c2r(solver->pm, workspace);
        
        fastpm_utils_readout(solver->pm, solver->p, workspace, get_position, ACC[d]);
        fastpm_info("d= %d ACC = %d\n", d, ACC[d]);
    }

    /* now we paint \Psi by the lagrangian position q */

    memset(Fk, 0, sizeof(Fk[0]) * pm_size(solver->pm));

    for(d = 0; d < 3; d ++) {
        fastpm_utils_paint(solver->pm, solver->p, NULL, workspace2, get_lagrangian_position, ACC[d]);
        fastpm_apply_za_hmc_force_transfer(solver->pm, workspace2, workspace, d);

        /* add HMC force component to to Fk */
        ptrdiff_t ind;
        for(ind = 0; ind < pm_size(solver->pm); ind ++) {
            /* Wang's magic factor of 2 in 1301.1348. 
             * We do not put it in in hmc_force_2lpt_transfer */
            Fk[ind] += 2 * workspace[ind]; 
        }
    }
    pm_free(solver->pm, workspace2);
    pm_free(solver->pm, workspace);
}


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();
    MPI_Comm comm = MPI_COMM_WORLD;
    printf("%td\n", sizeof(FastPMFloat));
    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);


    FastPMHMCZA * self = alloca(sizeof(FastPMHMCZA));
    fastpm_hmc_za_init(self, 128, 512, 0.292, 12.0, comm);

    FastPMFloat * sigma = pm_alloc(self->pm);
    FastPMFloat * Fk = pm_alloc(self->pm);

    FastPMFloat * rho_init_ktruth = pm_alloc(self->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(self->pm);

    FastPMFloat * rho_init_k = pm_alloc(self->pm);
    FastPMFloat * rho_final_x = pm_alloc(self->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5.0e6, /* FIXME: sigma8 ~ 0.8 */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };

    fastpm_utils_fill_deltak(self->pm, rho_init_ktruth, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);

    fastpm_hmc_za_evolve(self, rho_init_ktruth, rho_final_xtruth);

    fastpm_utils_dump(self->pm, "rho_init_ktruth.raw", rho_init_ktruth);
    fastpm_utils_dump(self->pm, "rho_final_xtruth.raw", rho_final_xtruth);

    ptrdiff_t mode = 2 * pm_ravel_o_index(self->pm, (ptrdiff_t[]) {2, 2, 2});

    double delta = 0.1;
    ptrdiff_t i;
    for(i = 0; i < pm_size(self->pm); i ++) {
        sigma[i] = 1.0;
    }
    fastpm_utils_fill_deltak(self->pm, rho_init_k, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += delta;

    fastpm_hmc_za_evolve(self, rho_init_k, rho_final_x);
    double chisq1 = fastpm_hmc_za_chisq(self, rho_final_xtruth, rho_final_x, sigma);

    fastpm_utils_fill_deltak(self->pm, rho_init_k, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += 1.1 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k, rho_final_x);
    double chisq2 = fastpm_hmc_za_chisq(self, rho_final_xtruth, rho_final_x, sigma);

    fastpm_utils_fill_deltak(self->pm, rho_init_k, 301, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += 1.05 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k, rho_final_x);
    fastpm_hmc_za_force(self, rho_final_xtruth, rho_final_x, sigma, Fk);

    fastpm_utils_dump(self->pm, "Fk.raw", Fk);
    double analytic = Fk[mode];
    double numerical = (chisq2 - chisq1) / (0.1 * delta);
    fastpm_info("analytic = %g, numerical = %g, rat = %g\n", analytic, numerical, numerical / analytic);

    pm_free(self->pm, rho_final_x);
    pm_free(self->pm, rho_init_k);
    pm_free(self->pm, rho_final_xtruth);
    pm_free(self->pm, rho_init_ktruth);
    pm_free(self->pm, Fk);
    pm_free(self->pm, sigma);

    fastpm_hmc_za_destroy(self);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

