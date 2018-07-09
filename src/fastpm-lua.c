#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <getopt.h>

#include <fastpm/logging.h>
#include "lua-config.h"

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

void usage()
{
        fastpm_raise(-1, 
"Usage: fastpm-lua parameterfile [...] \n"
"\n"
"if main function is defined in the parameter file, execute it.\n"
);
}
int main(int argc, char * argv[]) {
    char * confstr;
    fastpm_set_msg_handler(_non_mpi_msg_handler, (MPI_Comm) 0, NULL);

    char opt;
    extern int optind;
    extern char * optarg;
    int dump_schema = 0;
    while ((opt = getopt(argc, argv, "Hh?")) != -1) {
        switch(opt) {
            case 'H':
                dump_schema = 1;
            break;
            case 'h':
            case '?':
            default:
                usage();
            break;
        }
    }
    char * filename = argv[optind];

    char * error;
    if(dump_schema) {
        confstr = lua_config_parse("_help", filename, argc - 1, argv + 1, &error);
        if(confstr == NULL) {
            fastpm_raise(-1, "%s\n", error);
        } else {
            fastpm_info("Supported Parameters are: \n%s\n", confstr);
        }
        return 0;
    } else {
        confstr = lua_config_parse("_parse_runmain", filename, argc - 1, argv + 1, &error);
        if(confstr == NULL) {
            fastpm_raise(-1, "%s\n", error);
        }
        fastpm_info("Compiled parameters are: \n%s\n", confstr);
        LuaConfig * config;
        config = lua_config_new(confstr);
        lua_config_free(config);
        free(confstr);
    }
    return 0;
}

