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

/* Helper function s */
static int luaL_eval(lua_State * L, const char * string)
{
    /* Evaluate string based on the global namespace at stack top */
    char * s = alloca(strlen(string) + 20);
    sprintf(s, "return %s", string);
    luaL_loadstring(L, s);
    lua_pushvalue(L, -2);
    lua_setupvalue(L, -2, 1);
    return lua_pcall(L, 0, 1, 0);
}

static char * _strdup(const char * str)
{
    /* C99 */
    char * ret = malloc(strlen(str) + 1);
    strcpy(ret, str);
    return ret;
}

struct enum_entry {
    char * str;
    int value;
};

static int parse_enum(const char * str, struct enum_entry * enum_table)
{
    struct enum_entry * p;
    for(p = enum_table; p->str; p ++) {
        if(!strcmp(str, p->str)) {
            return p->value;
        }
    }
    int n = 10;
    for(p = enum_table; p->str; p ++) {
        n += strlen(p->str) + 10;
    }
    char * options = malloc(n);
    options[0] = 0;
    for(p = enum_table; p->str; p ++) {
        if(p != enum_table)
            strcat(options, ", ");
        strcat(options, "`");
        strcat(options, p->str);
        strcat(options, "`");
    }

    fastpm_raise(9999, "value `%s` is not recognized. Options are %s \n",
        str, options);
    return 0;
}

/* evals the expression with -1 as env */
static char * read_string(lua_State * L, const char * name) 
{
    luaL_eval(L, name);
    const char * val = lua_tostring(L, -1);
    lua_pop(L, 1);
    if(val) return _strdup(val);
    return NULL;
}
static int read_enum(lua_State * L, const char * name, struct enum_entry * enumtype) {
    luaL_eval(L, name);
    const char * val = lua_tostring(L, -1);
    lua_pop(L, 1);
    return parse_enum(val, enumtype);
}
static int read_boolean(lua_State * L, const char * name) {
    luaL_eval(L, name);
    int val = lua_toboolean(L, -1);
    lua_pop(L, 1);
    return val;
}
static int read_integer(lua_State * L, const char * name) {
    luaL_eval(L, name);
    int val = lua_tointeger(L, -1);
    lua_pop(L, 1);
    return val;
}
static int * read_array_integer(lua_State * L, const char * name, int * len) {
    luaL_eval(L, name);
    const int n = luaL_len(L, -1);
    int * array = (int*) malloc(sizeof(int) * n);
    int i;
    for(i = 1; i <= n; ++i) {
        lua_pushinteger(L, i);
        lua_gettable(L, -2);
        int x = lua_tointeger(L, -1);
        lua_pop(L,1);
        array[i-1] = x;
    }
    lua_pop(L, 1);
    *len = n;
    return array;
}
static double read_number(lua_State * L, const char * name) {
    luaL_eval(L, name);
    double val = lua_tonumber(L, -1);
    lua_pop(L, 1);
    return val;
}
static double * read_array_number(lua_State * L, const char * name, int * len) {
    luaL_eval(L, name);
    const int n = luaL_len(L, -1);
    double * array = (double*) malloc(sizeof(double) * n);
    int i;
    for(i = 1; i <= n; ++i) {
        lua_pushinteger(L, i);
        lua_gettable(L, -2);
        double x = lua_tonumber(L, -1);
        lua_pop(L,1);
        array[i-1] = x;
    }
    lua_pop(L, 1);
    *len = n;
    return array;
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

    lua_getglobal(L, "_runmain");
    char * real = _strdup(filename); //realpath(filename, NULL);
    lua_pushstring(L, real);
    lua_pushboolean(L, runmain);

    int i;

    for(i = 0; i < argc; i ++) {
        lua_pushstring(L, argv[i]);
    }
    free(real);
    if(lua_pcall(L, 2 + argc, 1, 0)) {
        fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
    }
    confstr = _strdup(lua_tostring(L, -1));
    lua_pop(L, 1);

    lua_close(L);
    return confstr;
}

#include "parameters.h"

void
loads_param(char * confstr, Parameters * param)
{
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if(luaL_eval(L, confstr)) {
        /* This shall never happen unless the dump library is brokean */
        fastpm_raise(-1, "Error: %s\n", lua_tostring(L, -1));
    }

    param->string = _strdup(confstr);
    param->nc = read_integer(L, "nc");
    param->boxsize = read_number(L, "boxsize");

    param->time_step = read_array_number(L, "time_step", &param->n_time_step);

    param->aout = read_array_number(L, "output_redshifts", &param->n_aout);
    /* convert z to a */
    int i;
    for(i = 0; i < param->n_aout; i ++) {
        param->aout[i] = 1 / (param->aout[i] + 1);
    }
    param->omega_m = read_number(L, "omega_m");
    param->h = read_number(L, "h");

    param->pm_nc_factor = read_array_integer(L, "pm_nc_factor", &param->n_pm_nc_factor);
    param->change_pm = read_array_number(L, "change_pm", &param->n_change_pm);
    param->np_alloc_factor = read_number(L, "np_alloc_factor");

    // File Names and optional parameters realated
    param->read_grafic= read_string(L, "read_grafic");
    param->read_noisek = read_string(L, "read_noisek");
    param->read_noise = read_string(L, "read_noise");
    param->read_runpbic = read_string(L, "read_runpbic");

    param->read_powerspectrum = read_string(L, "read_powerspectrum");
    param->sigma8 = read_number(L, "sigma8");
    param->random_seed = read_integer(L, "random_seed");
    param->inverted_ic = read_boolean(L, "inverted_ic");
    param->remove_cosmic_variance = read_boolean(L, "remove_cosmic_variance");

    param->write_powerspectrum = read_string(L, "write_powerspectrum");

    param->write_snapshot = read_string(L, "write_snapshot");
    param->write_runpb_snapshot = read_string(L, "write_runpb_snapshot");
    param->write_noisek = read_string(L, "write_noisek");
    param->write_noise = read_string(L, "write_noise");

    {
    struct enum_entry table[] = {
        {"cola", FORCE_MODE_COLA},
        {"pm", FORCE_MODE_PM},
        {"zola", FORCE_MODE_ZOLA},
        {NULL, -1},
    };

    param->force_mode = read_enum(L, "force_mode", table);
    }

    param->cola_stdda = read_boolean(L, "cola_stdda");
    {
    struct enum_entry table[] = {
        {"linear", MODEL_LINEAR},
        {"2lpt", MODEL_2LPT},
        {"za", MODEL_ZA},
        {"pm", MODEL_PM},
        {"none", MODEL_NONE},
        {NULL, -1},
    };
    param->enforce_broadband_mode = read_enum(L, "enforce_broadband_mode", table);
    }

    param->enforce_broadband_kmax = read_integer(L, "enforce_broadband_kmax");

    param->use_dx1_only = read_boolean(L, "za");

    {
    struct enum_entry table[] = {
        {"eastwood", KERNEL_EASTWOOD},
        {"3_4", KERNEL_3_4},
        {"5_4", KERNEL_5_4},
        {"3_2", KERNEL_3_2},
        {NULL, -1},
    };
    param->kernel_type = read_enum(L, "kernel_type", table);
    }
    lua_close(L);
}

int read_parameters(char * filename, Parameters * param, int argc, char ** argv, MPI_Comm comm)
{
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    /* read the configuration file */
    char * confstr;
    int confstr_len;

    /* run the parameter file on root rank.
     * other ranks use the serialized string to avoid duplicated
     * error reports */
    if(ThisTask == 0) {
        confstr = run_paramfile(filename, 0, argc, argv);
        confstr_len = strlen(confstr) + 1;
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    } else {
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        confstr = malloc(confstr_len);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    }

    fastpm_info("Configuration %s\n", confstr);

    loads_param(confstr, param);

    free(confstr);
    return 0;
}

