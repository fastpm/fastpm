#include "pmpfft.h"
#include "pm2lpt.h"
#include "vpm.h"
#include "walltime.h"

typedef double (*fastpm_pkfunc)(double k, void * data);

#ifdef __cplusplus
extern "C" {
#endif
int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm);

void
fastpm_init_pm(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void 
fastpm_evolve_2lpt(PM * pm, PMStore * pdata, 
        double a, double omega_m, 
        real_t * deltak_0);

void 
fastpm_derivative_2lpt(PM * pm, 
        PMStore * p, /* Current position (x) saved in -> x */
        real_t * rhopx, /* rho_p (the data) in x-space */
        real_t * Fk     /* (out) hmc force in fourier space */
        );

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

#ifdef __cplusplus
}
#endif
