#include <mpi.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#include "pmpfft.h"

void fastpm_utils_init_randtable();
FastPMMemory GMEM;

void libfastpm_init()
{
    pm_module_init();
    fastpm_utils_init_randtable();
    fastpm_set_msg_handler(fastpm_void_msg_handler, MPI_COMM_WORLD, NULL);
    GMEM.alignment = 1024 * 4;
    fastpm_memory_init(&GMEM, 0, 1);
}

void libfastpm_cleanup()
{
    fastpm_memory_destroy(&GMEM);
    pm_module_cleanup();
}

void libfastpm_set_memory_bound(size_t size, int allow_unordered)
{
    fastpm_memory_destroy(&GMEM);
    fastpm_memory_init(&GMEM, size, allow_unordered);
}

FastPMMemory * _libfastpm_get_gmem()
{
    return &GMEM;
}
