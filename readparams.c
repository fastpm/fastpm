#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "parameters.h"
#include "msg.h"
#include "fastpm-preface.h"


static void 
parse_conf(char * confstr, Parameters * param, lua_State * L);

#define DEF_READ(my_typename, c_typename, lua_typename, transform) \
c_typename read_## my_typename ## _opt (lua_State * L, const char * name, c_typename defvalue) { \
    lua_getglobal(L, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        return defvalue; \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1)); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename read_## my_typename (lua_State * L, const char * name) { \
    lua_getglobal(L, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        msg_abort(1030, "Error: Parameter %s not found in the parameter file\n", name); \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1)); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename * read_array_ ## my_typename(lua_State* L, const char * name, int *len) \
{ \
    lua_getglobal(L, name); \
    if(!lua_istable(L, -1)) { \
        msg_abort(1031, "Error: Parameter %s not found or not an array in the parameter file\n", name); \
    } \
    const int n = luaL_len(L, -1);     \
    c_typename * array = (c_typename *) malloc(sizeof(c_typename) * n); \
    int i; \
    for(i = 1; i <= n; ++i) { \
        lua_pushinteger(L, i); \
        lua_gettable(L, -2); \
        c_typename x = transform(lua_to ## lua_typename(L, -1)); \
        lua_pop(L,1); \
        array[i-1] = x; \
    } \
    lua_pop(L, 1); \
    *len = n; \
    return array; \
}

#define DEF_READ2(my_typename, c_typename, lua_typename, transform, argtype, argname) \
c_typename read_## my_typename ## _opt (lua_State * L, const char * name, c_typename defvalue, argtype argname){ \
    lua_getglobal(L, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        return defvalue; \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1), argname); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename read_## my_typename (lua_State * L, const char * name, argtype argname) { \
    lua_getglobal(L, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        msg_abort(1030, "Error: Parameter %s not found in the parameter file\n", name); \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1), argname); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename * read_array_ ## my_typename(lua_State* L, const char * name, int *len, argtype argname) \
{ \
    lua_getglobal(L, name); \
    if(!lua_istable(L, -1)) { \
        msg_abort(1031, "Error: Parameter %s not found or not an array in the parameter file\n", name); \
    } \
    int n = luaL_len(L, -1);     \
    c_typename * array = (c_typename *) malloc(sizeof(c_typename) * n); \
    int i; \
    for(i = 1; i <= n; ++i) { \
        lua_pushinteger(L, i); \
        lua_gettable(L, -2); \
        c_typename x = transform(lua_to ## lua_typename(L, -1), argname); \
        lua_pop(L,1); \
        array[i-1] = x; \
    } \
    lua_pop(L, 1); \
    *len = n; \
    return array; \
}

#define PLAIN(a) (a)
static char * _strdup(const char * str) {
    char * ret = malloc(strlen(str) + 1);
    strcpy(ret, str);
    return ret;
}
struct enum_entry {
    char * str;
    int value;
};
static int parse_enum(const char * str, struct enum_entry * enum_table) {
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
     
    msg_abort(9999, "value `%s` is not recognized. Options are %s \n",
        str, options);
    return 0;
}
DEF_READ(boolean, int, boolean, PLAIN)
DEF_READ(integer, int, integer, PLAIN)
DEF_READ(number, double, number, PLAIN)
DEF_READ(string, char *, string, _strdup)
DEF_READ2(enum, int, string, parse_enum, struct enum_entry *, enum_table)

int read_parameters(char * filename, Parameters * param, MPI_Comm comm)
{
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    /* Load the preface */
    char * preface = alloca(fastpm_preface_lua_len + 1);
    memcpy(preface, fastpm_preface_lua, fastpm_preface_lua_len);
    preface[fastpm_preface_lua_len] = 0;

    if(luaL_loadstring(L, preface) || lua_pcall(L, 0, 0, 0)) {
        msg_abort(-1, "%s\n", lua_tostring(L, -1));
        return -1;
    }

    /* read the configuration file */
    char * confstr;
    int confstr_len;
    if(ThisTask == 0) {
        lua_getglobal(L, "read_file");
        lua_pushstring(L, filename);
        if(lua_pcall(L, 1, 1, 0)) {
            msg_abort(-1, "%s\n", lua_tostring(L, -1));
        }
        confstr = strdup(lua_tostring(L, -1));
        confstr_len = strlen(confstr) + 1;
        lua_pop(L, 1);

        /* 0-th rank parse first to report possible errors */
        parse_conf(confstr, param, L);

        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    } else {
        /* other ranks parse after the collective call to avoid duplicated
         * error reports */
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        confstr = malloc(confstr_len);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);

        parse_conf(confstr, param, L);
    }

    free(confstr);
    lua_close(L);
    return 0;
}

static void 
parse_conf(char * confstr, Parameters * param, lua_State * L) 
{
    if(luaL_loadstring(L, confstr) || lua_pcall(L, 0, 0, 0)) {
        msg_abort(-1, "%s\n", lua_tostring(L, -1));
    }

    memset(param, 0, sizeof(*param));
    param->nc = read_integer(L, "nc");
    param->boxsize = read_number(L, "boxsize");

    double a_init = read_number_opt(L, "a_init", -1);
    if(a_init != -1) {
        double a_final = read_number(L, "a_final");
        int ntimestep = read_integer(L, "ntimestep");
        param->time_step = malloc(sizeof(double) * (ntimestep + 1));
        param->n_time_step = ntimestep + 1;
        int i = 0;
        for(i = 0; i <= ntimestep ; i ++) {
            param->time_step[i] = (a_final - a_init) / ntimestep * i + a_init;
        }
        param->time_step[ntimestep] = a_final;
    } else {
        param->time_step = read_array_number(L, "time_step", &param->n_time_step);
    }

    param->zout = read_array_number(L, "output_redshifts", &param->n_zout);

    param->random_seed = read_integer(L, "random_seed");
    param->nrealization = read_integer(L, "nrealization");

    param->omega_m = read_number(L, "omega_m");
    param->h = read_number(L, "h");
    param->sigma8 = read_number(L, "sigma8");

    param->pm_nc_factor = read_array_integer(L, "pm_nc_factor", &param->n_pm_nc_factor);
    param->change_pm = read_array_number(L, "change_pm", &param->n_change_pm);
    if(param->n_change_pm != param->n_pm_nc_factor) {
        msg_abort(1001, "Error: length of change_pm and pm_nc_factor mismatch.");
    }
    param->np_alloc_factor = read_number(L, "np_alloc_factor");
    param->loglevel = read_integer(L, "loglevel");

    // File Names and optional parameters realated
    param->readic_filename = read_string_opt(L, "readic", NULL);
    param->readnoise_filename = read_string_opt(L, "readnoise", NULL);

    if(!param->readic_filename) {
        param->power_spectrum_filename = read_string(L, "powerspectrum");
    }

    param->measure_power_spectrum_filename = read_string_opt(L, "measure_power", NULL);

    param->snapshot_filename = read_string_opt(L, "snapshot", NULL);

    struct enum_entry table[] = {
        {"a", TIME_STEP_A},
        {"loga", TIME_STEP_LOGA},
        {"growth", TIME_STEP_GROWTH},
        {"cola", FORCE_MODE_COLA},
        {"cola1", FORCE_MODE_COLA1},
        {"za", FORCE_MODE_ZA},
        {"2lpt", FORCE_MODE_2LPT},
        {"pm", FORCE_MODE_PM},
        {NULL, -1},
    };
    struct enum_entry mond_table[] = {
        {"none", PM_MOND_NONE},
        {"simple", PM_MOND_SIMPLE},
        {"nbc", PM_MOND_NBC},
        {NULL, -1},
    };
    param->pm_mond_mode = read_enum_opt(L, "pm_mond_mode", PM_MOND_NONE, mond_table);
    if(param->pm_mond_mode > PM_MOND_NONE) {
        param->pm_mond_parameters = read_array_number(L, "pm_mond_parameters", &param->n_pm_mond_parameters);
    }
    param->force_mode = read_enum(L, "force_mode", table);
    param->enforce_broadband = read_boolean(L, "enforce_broadband");
    param->diff_order = read_integer(L, "diff_order");
    param->poisson_order = read_integer_opt(L, "poisson_order", 1);

    param->cola_stdda = read_boolean_opt(L, "cola_stdda", param->cola_stdda);

    if(param->force_mode == FORCE_MODE_2LPT ||
       param->force_mode == FORCE_MODE_ZA) {
        if(param->n_time_step != 2) 
            msg_abort(-1, "only one step is supported in 2LPT and ZA mode\n");
    }
}

