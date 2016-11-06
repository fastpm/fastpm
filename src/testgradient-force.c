#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

static double
chi2(PM * pm, FastPMStore * p)
{
    double sum = 0;
    ptrdiff_t i = 0;
    for(i = 0; i < p->np; i ++) {
        sum += p->acc[i][0] * p->acc[i][0];
        sum += p->acc[i][1] * p->acc[i][1];
        sum += p->acc[i][2] * p->acc[i][2];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));

    return sum * 0.5;
}

static void
chi2_gradient(PM * pm, FastPMStore * p, FastPMStore * out)
{
    fastpm_store_copy(p, out);
}

static double
objective(FastPMSolver * solver,  FastPMStore * p)
{
    FastPMFloat * delta_k = pm_alloc(solver->basepm);

    fastpm_gravity_calculate(solver->gravity, solver->basepm, p, delta_k);

    double s = chi2(solver->basepm, p);

    pm_free(solver->basepm, delta_k);
    return s;
}

static void
objective_gradient(FastPMSolver * solver, FastPMStore * p, FastPMStore * grad_pos)
{
    FastPMStore tmp[1];
    fastpm_store_init(tmp);
    fastpm_store_alloc(tmp, p->np, p->attributes);

    FastPMFloat * delta_k = pm_alloc(solver->basepm);

    fastpm_gravity_calculate(solver->gravity, solver->basepm, p, delta_k);
    pm_free(solver->basepm, delta_k);
    chi2_gradient(solver->basepm, p, tmp);

    fastpm_gravity_calculate_gradient(solver->gravity, solver->basepm,
            tmp, p, grad_pos);

    fastpm_store_destroy(tmp);
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
        .PAINTER_TYPE = FASTPM_PAINTER_QUAD,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
    };

    FastPMSolver solver[1];
    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);
    FastPMStore gradient[1];

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 50000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };

    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    fastpm_solver_setup_ic(solver, rho_init_ktruth);
    fastpm_solver_evolve(solver, (double[]){0.1, }, 1);
    {
        //solver->p->x[0][0] = 1.4;
        //solver->p->x[0][1] = 0;
        //solver->p->x[0][2] = 0;
    }
    fastpm_store_init(gradient);
    fastpm_store_alloc(gradient, solver->p->np, solver->p->attributes);
    fastpm_store_copy(solver->p, gradient);


    objective_gradient(solver, solver->p, gradient);

    double obj1 = objective(solver, solver->p);

    int i, d;
    for(i = 0; i < 1; i ++) {
        fastpm_info("x[%d] = {%g %g %g }\n", i, 
            solver->p->x[i][0],
            solver->p->x[i][1],
            solver->p->x[i][2]);

        for(d = 0; d <= 2; d++) {
            double old = solver->p->x[i][d];

            solver->p->x[i][d] += 1e-5;

            double obj2 = objective(solver, solver->p);

            solver->p->x[i][d] = old;

            fastpm_info("objective1 = %.18e\n", obj1);
            fastpm_info("objective2 = %.18e\n", obj2);

            fastpm_info("gradient ad  = %g\n", gradient->x[i][d]);
            fastpm_info("gradient num = %g\n", (obj2 - obj1) / 1e-5);
            fastpm_info("gradient ratio  = %g\n", gradient->x[i][d] / ((obj2 - obj1) / 1e-5));
        }

    }
    fastpm_store_destroy(gradient);

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

