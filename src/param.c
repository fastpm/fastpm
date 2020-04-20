#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>

#include <getopt.h>

#include "lua-config.h"
#include "param.h"

static char * _strdup(const char * s)
{
    char * rt = malloc(strlen(s) + 1);
    strcpy(rt, s);
    return rt;
}
extern char * 
lua_config_parse(char * entrypoint, char * filename, int argc, char ** argv, char ** error);

CLIParameters *
parse_cli_args(int argc, char ** argv)
{
    CLIParameters * prr = malloc(sizeof(prr[0]));

    int opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    prr->NprocY = 0;
    prr->Nwriters = 0;
    prr->MemoryPerRank = 0;
    prr->MaxThreads = -1;
    prr->RestartSnapshotPath = NULL;
    while ((opt = getopt(argc, argv, "h?T:y:fW:m:r:")) != -1) {
        switch(opt) {
            case 'r':
                prr->RestartSnapshotPath = _strdup(optarg);
            break;
            case 'y':
                prr->NprocY = atoi(optarg);
            break;
            case 'f':
                prr->UseFFTW = 1;
            break;
            case 'T':
                prr->MaxThreads = atoi(optarg);
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

    prr->argv = argv + optind;
    prr->argc = argc - optind;

    return prr;

usage:
    printf("Usage: fastpm [-T MaxThreads] [-W Nwriters] [-f] [-y NprocY] [-m MemoryBoundInMB] paramfile\n"
    "-T limit number of OMP threads\n"
    "-f Use FFTW / slab decomposition \n"
    "-m limit memory usage (die if exceeds this)\n"
    "-y Set the number of processes in the 2D mesh along the Y direction. \n"
    "-r Restart from a given snapshot.\n"
);
    free(prr);
    return NULL;
}

LUAParameters *
parse_config(char * filename, int argc, char ** argv, char ** error)
{
    *error = NULL;

    LUAParameters * prr = malloc(sizeof(prr[0]));
    char * confstr = lua_config_parse("_parse", filename, argc, argv, error);
    if(confstr == NULL) {
        free(prr);
        return NULL;
    }

    prr->config = lua_config_new(confstr);

    if(lua_config_error(prr->config)) {
        const char * err = lua_config_error(prr->config);
        *error = _strdup(err);
        free(prr);
        return NULL;
    }

    prr->string = _strdup(confstr);
    free(confstr);
    return prr;
}

void
free_lua_parameters(LUAParameters * param)
{
    lua_config_free(param->config);
    free(param->string);
    free(param);
}

void
free_cli_parameters(CLIParameters * param)
{
    free(param);
}

