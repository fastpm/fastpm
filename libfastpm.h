#include "pmpfft.h"
#include "pm2lpt.h"

typedef double (*fastpm_pkfunc)(double k, void * data);

int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm);

void 
fastpm_evolve_2lpt(PM * pm, PMStore * pdata, 
        double a, double omega_m, 
        real_t * deltak_0, real_t * deltak_1, MPI_Comm comm);

void 
fastpm_derivative_2lpt(PM * pm, 
        PMStore * p, /* Current position (x) saved in -> x */
        real_t * rhod_k, /* rhod in fourier space */
        real_t * Fk,     /* (out) hmc force in fourier space */
        MPI_Comm comm);

void 
fastpm_fill_deltak(PM * pm, real_t * deltak, int seed, fastpm_pkfunc pk, void * pkdata);

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */
