#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

static double
chi2(PM * pm, FastPMFloat * rho_k)
{
    PMKIter kiter;
    double sum = 0;
    for(pm_kiter_init(pm, &kiter);
        !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        double w = 2.0;
        if(kiter.iabs[2] == 0) w = 1.0;
        if(kiter.iabs[2] == pm_nmesh(pm)[2] / 2 + 1) w = 1.0;
        sum += w * 0.5 * rho_k[kiter.ind] * rho_k[kiter.ind];
        sum += w * 0.5 * rho_k[kiter.ind + 1] * rho_k[kiter.ind + 1];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    return sum;
}

static void
chi2_gradient(PM * pm, FastPMFloat * rho_k, FastPMFloat * out)
{
    pm_assign(pm, rho_k, out);
}

static double
objective(FastPMSolver * solver, FastPMFloat * rho_init_x)
{
    FastPMFloat * copy = pm_alloc(solver->basepm);
    FastPMFloat * rho_init_k = pm_alloc(solver->basepm);

    pm_assign(solver->basepm, rho_init_x, copy);
    pm_r2c(solver->basepm, copy, rho_init_k);

    double s = chi2(solver->basepm, rho_init_k);

    pm_free(solver->basepm, rho_init_k);
    pm_free(solver->basepm, copy);
    return s;
}

static void
objective_gradient(FastPMSolver * solver, FastPMFloat * rho_init_x, FastPMFloat * out)
{
    FastPMFloat * tmp = pm_alloc(solver->basepm);

    /* redo rho_init_k, since we do not record rho_init_x on a tape */
    FastPMFloat * rho_init_k = pm_alloc(solver->basepm);

    pm_assign(solver->basepm, rho_init_x, tmp);
    pm_r2c(solver->basepm, tmp, rho_init_k);


    chi2_gradient(solver->basepm, rho_init_k, tmp);

    pm_r2c_gradient(solver->basepm, tmp, out);

    pm_free(solver->basepm, rho_init_k);
    pm_free(solver->basepm, tmp);
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

    FastPMFloat * gradient = pm_alloc(solver->basepm);
    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    if(1) {
        struct fastpm_powerspec_eh_params eh = {
            .Norm = 50000.0, /* FIXME: this is not any particular sigma8. */
            .hubble_param = 0.7,
            .omegam = 0.260,
            .omegab = 0.044,
        };

        fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);
    }
    fastpm_apply_modify_mode_transfer(solver->basepm, rho_init_ktruth, rho_init_ktruth, (ptrdiff_t[]) {0, 0, 0, 0}, 0.0);

    FastPMFloat * pert = pm_alloc(solver->basepm);

    int i;
    for(i = 0; i < 10; i ++) {
        pm_assign(solver->basepm, rho_init_ktruth, pert);
        pm_c2r(solver->basepm, pert);

        objective_gradient(solver, pert, gradient);

        pm_assign(solver->basepm, rho_init_ktruth, pert);
        pm_c2r(solver->basepm, pert);

        double obj1 = objective(solver, pert);
        fastpm_info("pert[i] = %g\n", pert[i]);

        pm_assign(solver->basepm, rho_init_ktruth, pert);
        pm_c2r(solver->basepm, pert);
        double v = 1e-3 * pert[i];

        pert[i] += v;

        double obj2 = objective(solver, pert);
        fastpm_info("----- i = %d ----\n", i);
        fastpm_info("obj1 = %.18e\n", obj1);
        fastpm_info("obj2 = %.18e\n", obj2);

        fastpm_info("gradient ad = %.18e\n", gradient[i]);
        fastpm_info("gradient num = %.18e\n", (obj2 - obj1) / v);
        fastpm_info("ratio  = %.18e\n", 1 / ((obj2 - obj1) / v / gradient[i]));
    }
    pm_free(solver->basepm, pert);
    pm_free(solver->basepm, rho_init_ktruth);
    pm_free(solver->basepm, gradient);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

