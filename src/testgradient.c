#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

static double
chi2(PM * pm, FastPMFloat * rho_x)
{
    PMXIter xiter;
    double sum = 0;
    for(pm_xiter_init(pm, &xiter);
        !pm_xiter_stop(&xiter);
        pm_xiter_next(&xiter)) {
        sum += 0.5 * rho_x[xiter.ind] * rho_x[xiter.ind];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    return sum;
}

static void
chi2_gradient(PM * pm, FastPMFloat * rho_x, FastPMFloat * out)
{
    pm_assign(pm, rho_x, out);
}

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 32,
        .boxsize = 32.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
    };

    FastPMSolver solver[1];
    fastpm_solver_init(solver, config, comm);

    FastPMFloat * tmp = pm_alloc(solver->basepm);
    FastPMFloat * grad1 = pm_alloc(solver->basepm);
    FastPMFloat * grad2 = pm_alloc(solver->basepm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);
    FastPMFloat * rho_init_xtruth = pm_alloc(solver->basepm);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->basepm);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    pm_assign(solver->basepm, rho_init_ktruth, tmp);
    pm_r2c(solver->basepm, tmp, rho_init_xtruth);

    double objective = chi2(solver->basepm, rho_init_xtruth);

    fastpm_info("objective = %g\n", objective);

    chi2_gradient(solver->basepm, rho_init_ktruth, grad1);
    pm_r2c_gradient(solver->basepm, grad1, grad2);

    pm_free(solver->basepm, rho_final_xtruth);
    pm_free(solver->basepm, rho_final_ktruth);
    pm_free(solver->basepm, rho_init_xtruth);
    pm_free(solver->basepm, rho_init_ktruth);
    pm_free(solver->basepm, grad2);
    pm_free(solver->basepm, grad1);
    pm_free(solver->basepm, tmp);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

