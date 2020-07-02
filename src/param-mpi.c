#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>

#include <getopt.h>
#include <fastpm/libfastpm.h>

#include "lua-config.h"
#include "param.h"

LUAParameters *
parse_config_mpi(char * filename, int argc, char ** argv, char ** error, MPI_Comm comm)
{
    int ThisTask;
    int NTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    LUAParameters * prr = NULL;

    *error = NULL;

    if(ThisTask == 0) {
        prr = parse_config(filename, argc, argv, error);
    } else {
        prr = malloc(sizeof(prr[0]));
    }

    if(MPIU_Any(comm, prr == NULL)) {
        if(prr != NULL) {
            free(prr);
        }

        if(MPIU_Any(comm, *error != NULL)) {
            *error = MPIU_Bcast_string(comm, *error, 0, free);
        }
        return NULL;
    }

    MPI_Bcast(prr, sizeof(prr[0]), MPI_BYTE, 0, comm);

    prr->string = MPIU_Bcast_string(comm, prr->string, 0, free);
    if(ThisTask == 0) {
        lua_config_free(prr->config);
    }
    prr->config = lua_config_new(prr->string);
    return prr;
}

CLIParameters *
parse_cli_args_mpi(int argc, char ** argv, MPI_Comm comm)
{
    CLIParameters * prr = parse_cli_args(argc, argv);

    int Nwriters = prr->Nwriters;
    if(Nwriters == 0) {
        MPI_Comm_size(comm, &Nwriters);
        prr->Nwriters = Nwriters;
    }
    return prr;
}
