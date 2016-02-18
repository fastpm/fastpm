#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <fastpm/logging.h>
#include "parameters.h"

extern void 
loads_param(char * confstr, Parameters * param, lua_State * L);

extern char * 
run_paramfile(char * filename, lua_State * L);

static void _non_mpi_msg_handler(
            const enum FastPMLogLevel level,
            const enum FastPMLogType type,
            const int errcode,
             const char * message, 
            MPI_Comm comm,
            void * userdata) {
    fprintf(stdout, "%s", message);
    fflush(stdout);

    if(level == ERROR) abort();

}

int main(int argc, char * argv[]) {
    char * confstr;
    fastpm_set_msg_handler(_non_mpi_msg_handler, (MPI_Comm) 0, NULL);

    if(argc < 2) {
        fastpm_raise(-1, "Usage: fastpm-validate filename\n");
    }
    int i;
    for(i = 1; i < argc; i ++) {
        char * filename = argv[i];

        lua_State *L = luaL_newstate();

        luaL_openlibs(L);

        confstr = run_paramfile(filename, L);

        Parameters param;
        loads_param(confstr, &param, L);

        free(confstr);
        lua_close(L);
    }
}

