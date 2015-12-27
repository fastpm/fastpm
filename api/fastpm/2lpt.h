FASTPM_BEGIN_DECLS

typedef struct FastPM2LPTSolver {
    PM * pm;
    PMStore * p;
    int Ngrid;
    MPI_Comm comm;
} FastPM2LPTSolver;

int 
fastpm_2lpt_init(FastPM2LPTSolver * solver, int nc, double BoxSize, double alloc_factor, MPI_Comm comm);

void 
fastpm_2lpt_evolve(FastPM2LPTSolver * solver, FastPMFloat * delta_k_i, double a, double omega_m);

void 
fastpm_2lpt_hmc_force(FastPM2LPTSolver * solver, 
        FastPMFloat * rhopx, /* (in), rho_p (the data) in x-space */
        FastPMFloat * Fk     /* (out) hmc force in fourier space */
        );

FASTPM_END_DECLS
