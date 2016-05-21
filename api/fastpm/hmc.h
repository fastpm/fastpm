#ifndef _FASTPM_HMC_H_
#define _FASTPM_HMC_H_

FASTPM_BEGIN_DECLS

typedef struct {
    double OmegaM;
    double SmoothingLength;
    double KThreshold;
    int DeconvolveCIC;
    int IncludeRSD;
    int LPTOrder;

    struct {
        fastpm_fkfunc func;
        void * data;
    } TransferFunction;

    int Nmesh;
    int Ngrid;
    double BoxSize;
    //int testing;

    PM * pm;
    PMStore * p;

    FastPM2LPTSolver solver;
    FastPM pm_solver;
    FastPMFloat * delta_ic_k;
    FastPMFloat * rho_final_x;
    FastPMFloat * transfer_function;

} FastPMHMCZA;

void
fastpm_hmc_za_init(FastPMHMCZA * self, MPI_Comm comm);


void
fastpm_hmc_za_destroy(FastPMHMCZA * self);

void
fastpm_hmc_za_evolve(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic, /* IC in k-space*/
    int Nsteps
    );

double
fastpm_hmc_za_chisq(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * sigma_x /* sigma_x in x-space*/
    );

void 
fastpm_hmc_za_force(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * Fk    /* (out) hmc force in fourier space */
    );

void
fastpm_hmc_za_force_rhodk(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * rhodk    /* (out) rhodk in fourier space */
    );

void
fastpm_hmc_za_force_s1(
    FastPMHMCZA * self,
    FastPMFloat * rhodk,
    FastPMFloat * Fk1    /* (out) hmc force for s2 in fourier space */
    );

void
fastpm_hmc_za_force_s2(
    FastPMHMCZA * self,
    FastPMFloat * Fk1,
    FastPMFloat * Fk2    /* (out) hmc force for s2 in fourier space */
    );

FASTPM_END_DECLS
#endif
