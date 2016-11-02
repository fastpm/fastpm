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
objective(FastPMSolver * solver,  FastPMStore * p)
{
    FastPMFloat * rho_x = pm_alloc(solver->basepm);
    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->basepm, FASTPM_PAINTER_CIC, 2);

    memset(rho_x, 0, sizeof(rho_x[0]) * pm_allocsize(solver->basepm));

    fastpm_paint(painter, rho_x, p, NULL, 0);

    double s = chi2(solver->basepm, rho_x);

    pm_free(solver->basepm, rho_x);
    return s;
}

static void
objective_gradient(FastPMSolver * solver, FastPMStore * p, FastPMStore * out)
{
    FastPMFloat * tmp = pm_alloc(solver->basepm);
    FastPMFloat * rho_x = pm_alloc(solver->basepm);

    /* redo rho_init_x, since we do not record rho_init_x on a tape */
    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->basepm, FASTPM_PAINTER_CIC, 2);

    memset(rho_x, 0, sizeof(rho_x[0]) * pm_allocsize(solver->basepm));

    fastpm_paint(painter, rho_x, p, NULL, 0);

    chi2_gradient(solver->basepm, rho_x, tmp);
    fastpm_paint_gradient(painter, tmp, p, NULL, 0, out);

    pm_free(solver->basepm, rho_x);
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
    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 50000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);
    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    fastpm_solver_evolve(solver, (double[]){1.0, }, 1);

    FastPMStore gradient[1];

    fastpm_store_init(gradient);
    fastpm_store_alloc(gradient, solver->p->np, solver->p->attributes);

    double obj1 = objective(solver, solver->p);

    solver->p->x[0][0] += 1e-5;
    double obj2 = objective(solver, solver->p);

    fastpm_info("objective1 = %.18e\n", obj1);
    fastpm_info("objective2 = %.18e\n", obj2);

    objective_gradient(solver, solver->p, gradient);
    fastpm_info("gradient ad  = %g\n", gradient->x[0][0]);
    fastpm_info("gradient num = %g\n", (obj2 - obj1) / 1e-5);
    fastpm_store_destroy(gradient);
    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

