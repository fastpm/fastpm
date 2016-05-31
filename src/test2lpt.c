#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/hmc.h>

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();
    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);


    FastPMHMCZA * self = &(FastPMHMCZA) {
            .BoxSize = 512.,
            .Nmesh = 128,
            .Ngrid = 128,
            .OmegaM = 0.304,
            .KThreshold = 0.0,
            .DeconvolveCIC = 0,
            .LPTOrder = 1,
            .SmoothingLength = 8,
            .IncludeRSD = 0,
	    .AllocFactor = 2.0,
        };

    fastpm_hmc_za_init(self, comm);

    FastPMFloat * sigma = pm_alloc(self->pm);

    FastPMFloat * rho_init_ktruth = pm_alloc(self->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(self->pm);

    FastPMFloat * rho_init_k0 = pm_alloc(self->pm);
    FastPMFloat * rho_init_k = pm_alloc(self->pm);

    FastPMFloat * Fk1 = pm_alloc(self->pm);
    FastPMFloat * Fk2 = pm_alloc(self->pm);
    FastPMFloat * rhodk = pm_alloc(self->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5.0e6, /* FIXME: sigma8 ~ 0.8 */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    ptrdiff_t mode = 2 * pm_ravel_o_index(self->pm, (ptrdiff_t[]) {1, 1, 2}) + 0;
    int seed = 299;

    fastpm_ic_fill_gaussiank(self->pm, rho_init_k0, seed, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(self->pm, rho_init_k0, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    pm_assign(self->pm, self->rho_final_x, rho_final_xtruth);
//    double amplitude = 10000;
//    fastpm_utils_fill_deltak(self->pm, rho_init_k0, seed, (fastpm_fkfunc)fastpm_utils_powerspec_white, &amplitude, FASTPM_DELTAK_GADGET);
    memset(rho_init_ktruth, 0, pm_size(self->pm) * sizeof(FastPMFloat));
//    fastpm_utils_fill_deltak(self->pm, rho_init_ktruth, seed, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);

    fastpm_hmc_za_evolve(self, rho_init_ktruth, 5);

    pm_assign(self->pm, self->rho_final_x, rho_final_xtruth);

    fastpm_utils_dump(self->pm, "rho_init_ktruth.raw", rho_init_ktruth);
    fastpm_utils_dump(self->pm, "rho_final_xtruth.raw", rho_final_xtruth);

    double delta = 1e-2;
    ptrdiff_t i;
    for(i = 0; i < pm_size(self->pm); i ++) {
        sigma[i] = 1.0;
    }
    /* Analytic */
    pm_assign(self->pm, rho_init_k0, rho_init_k);
    rho_init_k[mode] += 1.0005 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k, 5);
    fastpm_hmc_za_force_rhodk(self, rho_final_xtruth, sigma, rhodk);
    fastpm_hmc_za_force_s1(self, rhodk, Fk1);
    fastpm_hmc_za_force_s2(self, Fk1, Fk2);
    fastpm_utils_dump(self->pm, "Fk1.raw", Fk1);
    fastpm_utils_dump(self->pm, "Fk2.raw", Fk2);
    fastpm_utils_dump(self->pm, "rho_final_x.raw", self->rho_final_x);

    /* Numeric */
    pm_assign(self->pm, rho_init_k0, rho_init_k);
    rho_init_k[mode] += 1.0 * delta;
    fastpm_hmc_za_evolve(self, rho_init_k, 5);
    double chisq1 = fastpm_hmc_za_chisq(self, rho_final_xtruth, sigma);
    fastpm_utils_dump(self->pm, "rho_final_x1.raw", self->rho_final_x);

    pm_assign(self->pm, rho_init_k0, rho_init_k);
    rho_init_k[mode] += 1.001 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k, 5);
    double chisq2 = fastpm_hmc_za_chisq(self, rho_final_xtruth, sigma);
    fastpm_utils_dump(self->pm, "rho_final_x2.raw", self->rho_final_x);

    double analytic1 = Fk1[mode];
    double analytic2 = Fk2[mode];
    double analytic = Fk1[mode] + Fk2[mode];
    double numerical = (chisq2 - chisq1) / (0.001 * delta);

    fastpm_info("analytic1 = %g, numerical = %g, rat = %g\n", analytic1, numerical, analytic1 / numerical);
    fastpm_info("analytic2 = %g, numerical = %g, rat = %g\n", analytic2, numerical, analytic2 / numerical);
    fastpm_info("analytic = %g, numerical = %g, rat = %g\n", analytic, numerical, analytic / numerical);
    fastpm_info("analytic1[0] = %g\n", Fk1[0]);
    fastpm_info("analytic2[0] = %g\n", Fk2[0]);

    pm_free(self->pm, rhodk);
    pm_free(self->pm, Fk2);
    pm_free(self->pm, Fk1);

    pm_free(self->pm, rho_init_k);
    pm_free(self->pm, rho_init_k0);
    pm_free(self->pm, rho_final_xtruth);
    pm_free(self->pm, rho_init_ktruth);
    pm_free(self->pm, sigma);

    fastpm_hmc_za_destroy(self);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

