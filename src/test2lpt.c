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


    FastPMHMCZA * self = alloca(sizeof(FastPMHMCZA));
    fastpm_hmc_za_init(self, 256, 128, 512., 0.304, 0, comm);

    self->sml = 2.;
    self->kth = 0;
    self->decic = 1;

    FastPMFloat * sigma = pm_alloc(self->pm);
    FastPMFloat * Fk = pm_alloc(self->pm);

    FastPMFloat * rho_init_ktruth = pm_alloc(self->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(self->pm);

    FastPMFloat * rho_init_k = pm_alloc(self->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 5.0e6, /* FIXME: sigma8 ~ 0.8 */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    ptrdiff_t mode = 2 * pm_ravel_o_index(self->pm, (ptrdiff_t[]) {0, 0, 1}) + 0;
    int seed = 2999;

    fastpm_utils_fill_deltak(self->pm, rho_init_ktruth, seed, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    memset(rho_init_ktruth, 0, pm_size(self->pm) * sizeof(FastPMFloat));

    fastpm_hmc_za_evolve(self, rho_init_ktruth);

    pm_assign(self->pm, self->rho_final_x, rho_final_xtruth);

    fastpm_utils_dump(self->pm, "rho_init_ktruth.raw", rho_init_ktruth);
    fastpm_utils_dump(self->pm, "rho_final_xtruth.raw", rho_final_xtruth);

    double delta = 1e-2;
    ptrdiff_t i;
    for(i = 0; i < pm_size(self->pm); i ++) {
        sigma[i] = 1.0;
    }
    fastpm_utils_fill_deltak(self->pm, rho_init_k, seed, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += 1.0 * delta;
    fastpm_hmc_za_evolve(self, rho_init_k);
    double chisq1 = fastpm_hmc_za_chisq(self, rho_final_xtruth, sigma);

    fastpm_utils_fill_deltak(self->pm, rho_init_k, seed, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += 1.001 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k);
    double chisq2 = fastpm_hmc_za_chisq(self, rho_final_xtruth, sigma);

    fastpm_utils_fill_deltak(self->pm, rho_init_k, seed, (fastpm_pkfunc)fastpm_utils_powerspec_eh, &eh, FASTPM_DELTAK_GADGET);
    rho_init_k[mode] += 1.0005 * delta;

    fastpm_hmc_za_evolve(self, rho_init_k);
    fastpm_hmc_za_force(self, rho_final_xtruth, sigma, Fk);

    fastpm_utils_dump(self->pm, "Fk.raw", Fk);
    double analytic = Fk[mode];
    double numerical = (chisq2 - chisq1) / (0.001 * delta);
    fastpm_info("analytic = %g, numerical = %g, rat = %g\n", analytic, numerical, numerical / analytic);
    fastpm_info("analytic[0] = %g\n", Fk[0]);

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

