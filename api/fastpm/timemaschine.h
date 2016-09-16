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
} FastPMTEStep;

// Protos

int fastpm_tevo_block_len(FastPMTEEntry *template);

FastPMTEStates *fastpm_tevo_generate_states(FastPMTEStates *states, int cycles, FastPMTEEntry *template, double *ts);

void fastpm_tevo_destroy_states(FastPMTEStates *states);

double fastpm_tevo_i2t(double time_steps[], int i, int N);

void fastpm_tevo_print_states(FastPMTEStates *states);

FASTPM_END_DECLS
