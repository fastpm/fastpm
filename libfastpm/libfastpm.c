#include <mpi.h>
#include <fastpm/libfastpm.h>

#include "pmpfft.h"

void fastpm_utils_init_randtable();

void libfastpm_init() 
{
    pm_module_init();    
    fastpm_utils_init_randtable();
}

void libfastpm_cleanup() 
{
    pm_module_cleanup();    
}
