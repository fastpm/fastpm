#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <fastpm/logging.h>

static char * _strdup(const char * str)
{
    /* C99 */
    char * ret = malloc(strlen(str) + 1);
    strcpy(ret, str);
    return ret;
}

char *
run_paramfile(char * filename, int runmain, int argc, char ** argv) 
{

    char * confstr;

    extern int lua_open_runtime(lua_State * L);

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if(lua_open_runtime(L)) {
        fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
    }

    lua_getglobal(L, "_main");
    char * real = filename; //realpath(filename, NULL);
    lua_pushstring(L, real);
    lua_pushboolean(L, runmain);

    int i;

    for(i = 0; i < argc; i ++) {
        lua_pushstring(L, argv[i]);
    }

    if(lua_pcall(L, 2 + argc, 1, 0)) {
        fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
    }
    confstr = _strdup(lua_tostring(L, -1));
    lua_pop(L, 1);

    lua_close(L);
    return confstr;
}

