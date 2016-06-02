FASTPM_BEGIN_DECLS

typedef struct FastPM2LPTSolver {
    FastPMSolverBase base;

    int nc;
    double boxsize;
    int USE_SHIFT;
    int USE_DX1_ONLY;
} FastPM2LPTSolver;

int 
fastpm_2lpt_init(FastPM2LPTSolver * solver, int nmesh, int nc, double BoxSize, double alloc_factor, MPI_Comm comm);

void
fastpm_2lpt_destroy(FastPM2LPTSolver * solver);

void 
fastpm_2lpt_evolve(FastPM2LPTSolver * solver, FastPMFloat * delta_k_i, double a, double omega_m);

FASTPM_END_DECLS
