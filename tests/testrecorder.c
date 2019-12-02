#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

typedef struct {
    FastPMSolver * solver;
    FastPMStore * tape;
    int maxsteps;
    FastPMStore bare[1];
} FastPMRecorder;

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps);

void fastpm_recorder_record(FastPMRecorder * recorder, FastPMTransition * transition, FastPMStore * p);

void fastpm_recorder_seek(FastPMRecorder * recorder, FastPMState * state, FastPMStore * out);

void fastpm_recorder_destroy(FastPMRecorder * recorder);

void fastpm_recorder_gradient_dot(FastPMRecorder * recorder, FastPMStore * vector, FastPMStore * gradient);

static void
record_transition(FastPMSolver * solver, FastPMTransitionEvent * event, FastPMRecorder * recorder)
{
    fastpm_recorder_record(recorder, event->transition, solver->p);
}

static int 
target_func(void * pdata, ptrdiff_t i, void * data) 
{
    FastPMStore * p = (FastPMStore *) pdata;
    PM * pm = (PM*) data;
    double pos[3];
    fastpm_store_get_lagrangian_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
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
        .Omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
        .SAVE_Q = 1,
    };

    FastPMSolver solver[1];
    FastPMRecorder recorder[1];

    fastpm_solver_init(solver, config, comm);

    fastpm_solver_add_event_handler(solver, FASTPM_EVENT_TRANSITION,
            FASTPM_EVENT_STAGE_AFTER,
            (FastPMEventHandlerFunction) record_transition,
            recorder);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->pm);
    FastPMFloat * rho_final_ktruth = pm_alloc(solver->pm);
    FastPMFloat * rho_final_xtruth = pm_alloc(solver->pm);

    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->pm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->pm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
    fastpm_solver_setup_lpt(solver, rho_init_ktruth);

    fastpm_recorder_init(recorder, solver, 100);

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));

    FastPMPainter painter[1];
    fastpm_painter_init(painter, solver->pm, config->PAINTER_TYPE, config->painter_support);

    pm_clear(solver->pm, rho_final_xtruth);

    fastpm_paint(painter, rho_final_xtruth, solver->p, NULL, 0);
    //fastpm_utils_dump(solver->pm, "fastpm_rho_final_xtruth.raw", rho_final_xtruth);

    fastpm_recorder_destroy(recorder);

    pm_free(solver->pm, rho_final_xtruth);
    pm_free(solver->pm, rho_final_ktruth);
    pm_free(solver->pm, rho_init_ktruth);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps)
{
    /* Allocate FastPMRecorder for up to maxsteps states per x, v, q, acc*/
    recorder->maxsteps = maxsteps;
    recorder->solver = solver;
    recorder->tape = calloc(maxsteps, sizeof(FastPMStore));
    fastpm_store_init(recorder->bare, solver->p->np_upper, COLUMN_Q | COLUMN_ID);
    fastpm_store_copy(solver->p, recorder->bare);
    fastpm_store_decompose(recorder->bare, target_func, solver->pm, solver->comm);

    fastpm_store_sort_by_id(recorder->bare);

    size_t np_total = fastpm_store_get_np_total(solver->p, solver->comm);
    int step;
    for(step = 0; step < maxsteps; step ++) {
        FastPMStore * p = &recorder->tape[step];
        fastpm_store_init_evenly(p, np_total,
            COLUMN_POS | COLUMN_VEL | COLUMN_ACC, 1.1, solver->comm);
        p->np = recorder->bare->np;
        /* steal a reference from bare. Need to nullify the pointer before destroy the stores */
        p->q = recorder->bare->q;
        p->id = recorder->bare->id;
    }
}

void fastpm_recorder_record(FastPMRecorder * recorder, FastPMTransition * transition, FastPMStore * p)
{
    /* according to the initial and final state of the transition, add / modify records in recorder.
     * initialize new FastPMStore objects with COLUMN_X | COLUMN_V | COLUMN_Q | COLUMN_ACC if necessary.
     *
     * Important: We want to sort the particle into a particular order before saving them.
     * */

    FastPMStore tmp[1];

    switch(transition->action) {
        case FASTPM_ACTION_DRIFT:
            fastpm_store_init(tmp, p->np_upper, COLUMN_POS | COLUMN_ID | COLUMN_Q);
        break;
        case FASTPM_ACTION_KICK:
            fastpm_store_init(tmp, p->np_upper, COLUMN_VEL | COLUMN_ID | COLUMN_Q);
        break;
        case FASTPM_ACTION_FORCE:
            fastpm_store_init(tmp, p->np_upper, COLUMN_ACC | COLUMN_ID | COLUMN_Q);
        break;
    }
    fastpm_store_copy(p, tmp);

    /* move particles by their initial position */
    fastpm_store_decompose(tmp, target_func, recorder->solver->pm, recorder->solver->comm);

    /* sort p locally by ID to ensure consistency in position */
    fastpm_store_sort_by_id(tmp);

    FastPMStore * dst = &recorder->tape[transition->i.f];

    /* fastpm_store_local_sort(p, key_by_id) */
    if(dst->np != tmp->np) {
        fastpm_raise(-1, "the domain decompostion by lagrangian coordinate is broken\n");
    }

    switch(transition->action) {
        case FASTPM_ACTION_DRIFT:
            memcpy(dst->x, tmp->x, sizeof(tmp->x[0]) * tmp->np);
            p->a_x = transition->a.f;
        break;
        case FASTPM_ACTION_KICK:
            memcpy(dst->v, tmp->v, sizeof(tmp->v[0]) * tmp->np);
            p->a_v = transition->a.f;
        break;
        case FASTPM_ACTION_FORCE:
            memcpy(dst->acc, tmp->acc, sizeof(tmp->acc[0]) * tmp->np);
            /* FIXME: a_acc is not tracked in FastPMStore */
        break;
    }
    fastpm_store_destroy(tmp);
}

void fastpm_recorder_seek(FastPMRecorder * recorder, FastPMState * state, FastPMStore * out)
{
    /* modify x, v, p pointers of out */
    out->x = recorder->tape[state->x].x;
    out->v = recorder->tape[state->v].v;
    out->acc = recorder->tape[state->force].acc;
}

void fastpm_recorder_destroy(FastPMRecorder * recorder)
{
    int step;
    for(step = recorder->maxsteps - 1; step >= 0; step --) {
        FastPMStore * p = &recorder->tape[step];
        p->q = NULL;
        p->id = NULL;
        fastpm_store_destroy(p);
    }
    fastpm_store_destroy(recorder->bare);
    free(recorder->tape);
}

void fastpm_recorder_gradient(FastPMRecorder * recorder, FastPMStore * vector, FastPMStore * gradient)
{
    /* backward differentiation;
     * vector: the gradient from xi^2 to final particle position.
     * gradient : the gradieint from xi^2 to initial position.
     */

}
