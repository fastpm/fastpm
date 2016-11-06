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

static double
objective(FastPMSolver * solver, FastPMFloat * rho_init_k)
{
    FastPMFloat * rho_init_x = pm_alloc(solver->basepm);

    pm_assign(solver->basepm, rho_init_k, rho_init_x);
    pm_c2r(solver->basepm, rho_init_x);

    double s = chi2(solver->basepm, rho_init_x);

    pm_free(solver->basepm, rho_init_x);
    return s;
}

static void
objective_gradient(FastPMSolver * solver, FastPMFloat * rho_init_k, FastPMFloat * out)
{
    FastPMFloat * tmp = pm_alloc(solver->basepm);
    FastPMFloat * rho_init_x = pm_alloc(solver->basepm);

    /* redo rho_init_x, since we do not record rho_init_x on a tape */
    pm_assign(solver->basepm, rho_init_k, rho_init_x);
    pm_c2r(solver->basepm, rho_init_x);

    chi2_gradient(solver->basepm, rho_init_x, tmp);

    pm_c2r_gradient(solver->basepm, tmp, out);

    pm_compress_gradient(solver->basepm, out, out);
    pm_free(solver->basepm, rho_init_x);
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
    FastPMFloat * rho_init_kperturb = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 50000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    fastpm_apply_modify_mode_transfer(solver->basepm, rho_init_ktruth, rho_init_ktruth, (ptrdiff_t[]) {0, 0, 0, 0}, 0.0);
    ptrdiff_t modes[][4] = {
            {1, 1, 1, 0,},
            {0, 1, 1, 0,},
            {0, 1, 1, 1,},
            {1, 0, 1, 0,},
            {1, 0, 0, 0,},
            {1, 0, 0, 1,},
            {-1, -1, -1, -1},
    };

    double obj = objective(solver, rho_init_ktruth);
    fastpm_info("objective = %g\n", obj);

    objective_gradient(solver, rho_init_ktruth, gradient);

    int i;
    for(i = 0; modes[i][0] >= 0; i ++) {
        double value = fastpm_apply_get_mode_transfer(solver->basepm, rho_init_ktruth, modes[i]);
        fastpm_info("mode = %td %td %td %td \n", modes[i][0], modes[i][1], modes[i][2], modes[i][3]);
        fastpm_apply_modify_mode_transfer(solver->basepm, rho_init_ktruth, rho_init_kperturb, modes[i], value + 1e-5);
        double value1 = fastpm_apply_get_mode_transfer(solver->basepm, rho_init_kperturb, modes[i]);

        double obj1 = objective(solver, rho_init_kperturb);
        fastpm_info("value = %g obj = %g\n", value, obj);
        fastpm_info("value1 = %g obj1 = %g\n", value1, obj1);

        double gnum = (obj1 - obj) / (value1 - value);
        double gany = fastpm_apply_get_mode_transfer(solver->basepm, gradient, modes[i]);
        fastpm_info("numerical gradient = %g\n", gnum);
        fastpm_info("analytic  gradient = %g\n", gany);
        fastpm_info("ratio = %g\n", gany / gnum);
    }


    pm_free(solver->basepm, rho_init_kperturb);
    pm_free(solver->basepm, rho_init_ktruth);
    pm_free(solver->basepm, gradient);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

