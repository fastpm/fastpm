FASTPM_BEGIN_DECLS

// Typedefs

typedef struct {
    int force, x, v;
} FastPMState;

typedef struct {
    FastPMState * table;
    int cycle_len;
    int cycles;
    double * timesteps;
} FastPMStates;

struct FastPMTransition {
    FastPMStates * states;
    int istart;
    int iend;
    FastPMState * start;
    FastPMState * end;
    enum FastPMAction action;
    struct {
        double i, f, r;
    } a;
    struct {
        int i, f, r;
    } i;
};

// Protos

FastPMStates * fastpm_tevo_generate_states(FastPMStates *states, int cycles, FastPMState * templ, double * ts);

void fastpm_tevo_destroy_states(FastPMStates *states);

void
fastpm_tevo_transition_init(FastPMTransition * transition, FastPMStates * states, int istart, int iend);

int
fastpm_tevo_transition_find_dual(FastPMTransition * transition, FastPMTransition * dual);

int
fastpm_tevo_transition_find_next(FastPMTransition * transition, FastPMTransition * next);

FASTPM_END_DECLS
