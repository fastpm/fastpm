#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <fastpm/logging.h>
#include "parameters.h"

extern void 
loads_param(char * confstr, Parameters * param);

extern char * 
run_paramfile(char * filename, int runmain, int argc, char ** argv);

static void _non_mpi_msg_handler(
            const enum FastPMLogLevel level,
            const enum FastPMLogType type,
            const int errcode,
             const char * message, 
            MPI_Comm comm,
            void * userdata) {
    fprintf(stdout, "%s", message);
    fflush(stdout);

    if(level == ERROR) exit(1);

}

int main(int argc, char * argv[]) {
    char * confstr;
    fastpm_set_msg_handler(_non_mpi_msg_handler, (MPI_Comm) 0, NULL);

    if(argc < 2) {
        fastpm_raise(-1, 
"Usage: fastpm-lua parameterfile [...] \n"
"\n"
"if main function is defined in the parameter file, execute it.\n"
);
    }
    char * filename = argv[1];

    confstr = run_paramfile(filename, 1, argc - 1, argv + 1);

    Parameters param;
    loads_param(confstr, &param);

    free(confstr);
}

