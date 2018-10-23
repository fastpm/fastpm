#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>

#include <getopt.h>

#include "lua-config.h"
#include "param.h"

extern char * 
lua_config_parse(char * entrypoint, char * filename, int argc, char ** argv, char ** error);

Parameters *
parse_args(int argc, char ** argv, char ** error)
{
    char * ParamFileName;
    Parameters * prr = malloc(sizeof(prr[0]));
    *error = NULL;
    int opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    ParamFileName = NULL;
    prr->NprocY = 0;
    prr->Nwriters = 0;
    prr->MemoryPerRank = 0;
    while ((opt = getopt(argc, argv, "h?y:fW:m:")) != -1) {
        switch(opt) {
            case 'y':
                prr->NprocY = atoi(optarg);
            break;
            case 'f':
                prr->UseFFTW = 1;
            break;
            case 'W':
                prr->Nwriters = atoi(optarg);
            break;
            case 'm':
                prr->MemoryPerRank = atoi(optarg);
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

    ParamFileName = (argv)[optind];
    argv += optind;
    argc -= optind; 

    char * confstr = lua_config_parse("_parse", ParamFileName, argc, argv, error);
    if(confstr == NULL) {
        free(prr);
        return NULL;
    }

    prr->config = lua_config_new(confstr);

    if(lua_config_error(prr->config)) {
        *error = strdup(lua_config_error(prr->config));
        free(prr);
        return NULL;
    }

    prr->string = strdup(confstr);
    free(confstr);

    return prr;

usage:
    printf("Usage: fastpm [-W Nwriters] [-f] [-y NprocY] [-m MemoryBoundInMB] paramfile\n"
    "-f Use FFTW \n"
    "-y Set the number of processes in the 2D mesh\n"
    "-n Throttle IO (bigfile only) \n"
);
    free(prr);
    return NULL;
}

void
free_parameters(Parameters * param)
{
    lua_config_free(param->config);
    free(param->string);
}

