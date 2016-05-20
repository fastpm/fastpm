#include <mpi.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#include "pmpfft.h"

void fastpm_utils_init_randtable();

void libfastpm_init()
{
    pm_module_init();
    fastpm_utils_init_randtable();
    fastpm_set_msg_handler(fastpm_void_msg_handler, MPI_COMM_WORLD, NULL);
}

void libfastpm_cleanup()
{
    pm_module_cleanup();
}
