typedef double (*fastpm_pkfunc)(double k, void * data);

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
fastpm_2lpt_paint(FastPM2LPTSolver * solver, FastPMFloat * delta_x, FastPMFloat * delta_k);
        
void 
fastpm_2lpt_hmc_force(FastPM2LPTSolver * solver, 
        FastPMFloat * rhopx, /* (in), rho_p (the data) in x-space */
        FastPMFloat * Fk     /* (out) hmc force in fourier space */
        );

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */

void 
fastpm_fill_deltak(PM * pm, FastPMFloat * deltak, 
        int seed, fastpm_pkfunc pk, void * pkdata);

void fastpm_2lpt_dump(FastPM2LPTSolver * solver, char * filename, FastPMFloat *data);
