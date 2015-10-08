#include "pmpfft.h"
#include "pm2lpt.h"

int 
fastpm_init(PMStore * p, int nc, double alloc_factor, MPI_Comm comm);

int fastpm_particle_to_mesh(PM * pm, PMStore * p);
