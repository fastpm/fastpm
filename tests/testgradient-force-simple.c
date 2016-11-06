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
objective(FastPMSolver * solver,  FastPMStore * p, FastPMFloat * rho_2)
{
    FastPMFloat * rho_x = pm_alloc(solver->basepm);

    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->basepm, FASTPM_PAINTER_QUAD, 2);

    FastPMPainter painter2[1];
    fastpm_painter_init(painter2, solver->basepm, FASTPM_PAINTER_QUAD, 2);

    fastpm_paint(painter, rho_x, p, NULL, 0);

//    fastpm_readout(painter2, rho_x, p, fastpm_store_get_lagrangian_position, PACK_ACC_X);
    fastpm_readout(painter2, rho_x, p, NULL, PACK_ACC_X);

    double s = chi2(solver->basepm, p);

    pm_free(solver->basepm, rho_x);
    return s;
}

static void
objective_gradient(FastPMSolver * solver, FastPMStore * p, FastPMStore * grad_pos, FastPMFloat * rho_2)
{
    ptrdiff_t i;
    for(i = 0; i < grad_pos->np; i ++) {
        grad_pos->x[i][0] = 0;
        grad_pos->x[i][1] = 0;
        grad_pos->x[i][2]= 0;
    }

    FastPMStore grad_acc[1];
    fastpm_store_init(grad_acc);
    fastpm_store_alloc(grad_acc, p->np, p->attributes);

    FastPMStore grad_pos1[1];
    fastpm_store_init(grad_pos1);
    fastpm_store_alloc(grad_pos1, p->np, p->attributes);

    FastPMFloat * rho_x = pm_alloc(solver->basepm);
    FastPMFloat * grad_field= pm_alloc(solver->basepm);

    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->basepm, FASTPM_PAINTER_QUAD, 2);

    FastPMPainter painter2[1];
    fastpm_painter_init(painter2, solver->basepm, FASTPM_PAINTER_QUAD, 2);

    fastpm_paint(painter, rho_x, p, NULL, 0);

//    fastpm_readout(painter2, rho_x, p, fastpm_store_get_lagrangian_position, PACK_ACC_X);
    fastpm_readout(painter2, rho_x, p, NULL, PACK_ACC_X);

    chi2_gradient(solver->basepm, p, grad_acc);

//    fastpm_readout_gradient(painter2, grad_acc, rho_x, p, fastpm_store_get_lagrangian_position, PACK_ACC_X, grad_field, grad_pos1);
    fastpm_readout_gradient(painter2, grad_acc, rho_x, p, NULL, PACK_ACC_X, grad_field, grad_pos1);
    fastpm_info("%d : %g %g %g\n", 0,
            grad_pos1->x[0][0],
            grad_pos1->x[0][1],
            grad_pos1->x[0][2]);
    for(i = 0; i < grad_pos->np; i ++) {
        grad_pos->x[i][0] += grad_pos1->x[i][0];
        grad_pos->x[i][1] += grad_pos1->x[i][1];
        grad_pos->x[i][2] += grad_pos1->x[i][2];
    }

    fastpm_paint_gradient(painter, grad_field, p, NULL, 0, grad_pos1, NULL);
    fastpm_info("rho_x[0] = %g\n", rho_x[0]);
    fastpm_info("grad[0] = %g\n", grad_field[0]);

    for(i = 0; i < grad_pos->np; i ++) {
        grad_pos->x[i][0] += grad_pos1->x[i][0];
        grad_pos->x[i][1] += grad_pos1->x[i][1];
        grad_pos->x[i][2] += grad_pos1->x[i][2];
    }
    fastpm_info("%d : %g %g %g\n", 4,
            grad_pos1->x[0][0],
            grad_pos1->x[0][1],
            grad_pos1->x[0][2]);

    pm_free(solver->basepm, grad_field);
    pm_free(solver->basepm, rho_x);

    fastpm_store_destroy(grad_pos1);
    fastpm_store_destroy(grad_acc);
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
        .PAINTER_TYPE = FASTPM_PAINTER_CIC,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
        .SAVE_Q = 1,
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
        //solver->p->x[0][0] = 3.2;
        //solver->p->x[0][1] = 3.2;
        //solver->p->x[0][2] = 3.2;
    }
    fastpm_store_init(gradient);
    fastpm_store_alloc(gradient, solver->p->np, solver->p->attributes);
    fastpm_store_copy(solver->p, gradient);

    pm_c2r(solver->basepm, rho_init_ktruth);

    objective_gradient(solver, solver->p, gradient, rho_init_ktruth);

    double obj1 = objective(solver, solver->p, rho_init_ktruth);

    int i, d;
    for(i = 0; i < 1; i ++) {
        fastpm_info("x[%d] = [%g %g %g]\n", i,
            solver->p->x[i][0],
            solver->p->x[i][1],
            solver->p->x[i][2]);

        for(d = 0; d <= 2; d++) {
            double old = solver->p->x[i][d];

            solver->p->x[i][d] += 1e-5;

            double obj2 = objective(solver, solver->p, rho_init_ktruth);

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

