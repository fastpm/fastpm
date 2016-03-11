#ifndef _FASTPM_HMC_H_
#define _FASTPM_HMC_H_

FASTPM_BEGIN_DECLS

typedef struct {
    FastPM2LPTSolver solver;
    double OmegaM;
    double sml;
    double kth;
    int decic;
    PM * pm;
} FastPMHMCZA;

void 
fastpm_hmc_za_init(FastPMHMCZA * self, 
    int nc,
    double boxsize,
    double OmegaM, 
    MPI_Comm comm);

void 
fastpm_hmc_za_destroy(FastPMHMCZA * self);

void 
fastpm_hmc_za_evolve(
    FastPMHMCZA * self,
    FastPMFloat * delta_ic, /* IC in k-space*/
    FastPMFloat * delta_final /* final in x-space*/
    );

double
fastpm_hmc_za_chisq(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * model_x, /* rhop in x-space*/
    FastPMFloat * sigma_x /* sigma_x in x-space*/
    );

void 
fastpm_hmc_za_force(
    FastPMHMCZA * self,
    FastPMFloat * data_x, /* rhop in x-space*/
    FastPMFloat * model_x, /* rhop in x-space*/
    FastPMFloat * sigma_x, /* sigma_x in x-space*/
    FastPMFloat * Fk    /* (out) hmc force in fourier space */
    );

FASTPM_END_DECLS
#endif
