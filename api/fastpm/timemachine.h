FASTPM_BEGIN_DECLS

// Typedefs

typedef struct {
    int a, x, v;
} FastPMTEEntry;

typedef struct {
    FastPMTEEntry *table;
    int cycle_len;
    int cycles;
    double *timesteps;
} FastPMTEStates;

typedef struct {
    int i, f, r;
    double a_i, a_f, a_r;
} FastPMTETransition;

// Protos

FastPMTEStates * fastpm_tevo_generate_states(FastPMTEStates *states, int cycles, FastPMTEEntry *template, double *ts);

void fastpm_tevo_destroy_states(FastPMTEStates *states);

void
fastpm_tevo_transition_init(FastPMTETransition * transition, FastPMTEStates * states, int i, int r, int f);
void fastpm_tevo_print_states(FastPMTEStates *states);

FASTPM_END_DECLS
