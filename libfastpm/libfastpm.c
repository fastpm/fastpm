#include <mpi.h>
#include <fastpm/libfastpm.h>

#include "pmpfft.h"

void libfastpm_init() 
{
    pm_module_init();    
}

void libfastpm_cleanup() 
{
    pm_module_cleanup();    
}
