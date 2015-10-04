#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <mpi.h>
#include <signal.h>
#include <getopt.h>

#include "parameters.h"
#include "readparams.h"

/* command-line arguments */
static char * ParamFileName;

static void 
parse_args(int argc, char ** argv, Parameters * prr);

int fastpm(Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    MPI_Init(&argc, &argv);

    Parameters prr;

    parse_args(argc, argv, &prr);

    read_parameters(ParamFileName, &prr);

    MPI_Comm comm = MPI_COMM_WORLD; 

    fastpm(&prr, comm);

    pfft_cleanup();
    MPI_Finalize();
}

static void 
parse_args(int argc, char ** argv, Parameters * prr) 
{
    char opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    ParamFileName = NULL;
    prr->NprocY = 0;    
    while ((opt = getopt(argc, argv, "h?y:f")) != -1) {
        switch(opt) {
            case 'y':
                prr->NprocY = atoi(optarg);
            break;
            case 'f':
                prr->UseFFTW = 1;
            break;
            case 'h':
            case '?':
            default:
                goto usage;
            break;
        }
    }
    if(optind >= argc) {
        goto usage;
    }

    ParamFileName = argv[optind];
    return;

usage:
    msg_printf(-1, "Usage: fastpm [-f] [-y NprocY] paramfile\n"
    "-f Use FFTW \n"
    "-y Set the number of processes in the 2D mesh\n"
);
    MPI_Finalize();
    exit(1);
}
