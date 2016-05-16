/*
 * Interface FastPMPowerSpectrum with fastpm.c
 *
 * This file will go away.
 *
 * */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/powerspectrum.h>
#include <fastpm/string.h>

#include "power.h"

static FastPMPowerSpectrum ps;

void power_init(const char filename[], const double sigma8, MPI_Comm comm)
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    char * content;
    if(myrank == 0) {
        content = fastpm_file_get_content(filename);
        int size = strlen(content);
        MPI_Bcast(&size, 1, MPI_INT, 0, comm);
        MPI_Bcast(content, size + 1, MPI_BYTE, 0, comm);
    } else {
        int size = 0;
        MPI_Bcast(&size, 1, MPI_INT, 0, comm);
        content = malloc(size + 1);
        MPI_Bcast(content, size + 1, MPI_BYTE, 0, comm);
    }
    if (0 != fastpm_powerspectrum_init_from_string(&ps, content)) {
        fastpm_raise(-1, "Failed to parse the powerspectrum\n");
    }
    free(content);

    fastpm_info("Found %d pairs of values in input spectrum table\n", ps.size);

    double sigma8_input= fastpm_powerspectrum_sigma(&ps, 8);
    printf("%g \n", sigma8_input);
    fastpm_info("Input power spectrum sigma8 %f\n", sigma8_input);

    if(sigma8 > 0) {
        fastpm_info("Expected power spectrum sigma8 %g; correction applied. \n", sigma8);
        fastpm_powerspectrum_scale(&ps, pow(sigma8 / sigma8_input, 2));
    }
}

double PowerSpecWithData(double k, void * data) {
    return fastpm_powerspectrum_eval(&ps, k);
}

double PowerSpec(const double k)
{
    return fastpm_powerspectrum_eval(&ps, k);
}
