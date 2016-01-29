#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <fastpm/logging.h>
#include "parameters.h"


static void 
loads(char * confstr, Parameters * param, lua_State * L);

#define DEF_READ(my_typename, c_typename, lua_typename, transform) \
c_typename read_## my_typename ## _opt (lua_State * L, const char * name, c_typename defvalue) { \
    lua_getfield(L, -1, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        lua_pop(L, 1); \
        return defvalue; \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1)); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename read_## my_typename (lua_State * L, const char * name) { \
    lua_getfield(L, -1, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        fastpm_raise(1030, "Error: Parameter %s not found in the parameter file\n", name); \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1)); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename * read_array_ ## my_typename(lua_State* L, const char * name, int *len) \
{ \
    lua_getfield(L, -1, name); \
    if(!lua_istable(L, -1)) { \
        fastpm_raise(1031, "Error: Parameter %s not found or not an array in the parameter file\n", name); \
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
    lua_getfield(L, -1, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        lua_pop(L, 1); \
        return defvalue; \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1), argname); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename read_## my_typename (lua_State * L, const char * name, argtype argname) { \
    lua_getfield(L, -1, name); \
    if (!lua_is ## lua_typename (L, -1)) { \
        fastpm_raise(1030, "Error: Parameter %s not found in the parameter file\n", name); \
    } \
    c_typename retvalue = transform(lua_to ## lua_typename(L, -1), argname); \
    lua_pop(L, 1); \
    return  retvalue; \
} \
c_typename * read_array_ ## my_typename(lua_State* L, const char * name, int *len, argtype argname) \
{ \
    lua_getfield(L, -1, name); \
    if(!lua_istable(L, -1)) { \
        fastpm_raise(1031, "Error: Parameter %s not found or not an array in the parameter file\n", name); \
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
     
    fastpm_raise(9999, "value `%s` is not recognized. Options are %s \n",
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
    extern char lua_preface_lua[];
    extern unsigned int lua_preface_lua_len;
    extern char lua_dump_lua[];
    extern unsigned int lua_dump_lua_len;

    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    /* read the configuration file */
    char * confstr;
    int confstr_len;

    if(ThisTask == 0) {

        if(luaL_loadbuffer(L, lua_dump_lua, lua_dump_lua_len, "dump") 
        || lua_pcall(L, 0, 1, 0)
        ) {
            fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
            return -1;
        }
        lua_setglobal(L, "dump");

        if(luaL_loadbuffer(L, lua_preface_lua, lua_preface_lua_len, "preface") 
        || lua_pcall(L, 0, 0, 0)) {
            fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
            return -1;
        }

        lua_getglobal(L, "parse_file");
        lua_pushstring(L, filename);
        if(lua_pcall(L, 1, 1, 0)) {
            fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
        }
        confstr = strdup(lua_tostring(L, -1));
        confstr_len = strlen(confstr) + 1;
        lua_pop(L, 1);
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    } else {
        /* other ranks use the serialized string to avoid duplicated
         * error reports */
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        confstr = malloc(confstr_len);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    }

    loads(confstr, param, L);

    free(confstr);
    lua_close(L);
    return 0;
}

static void 
loads(char * confstr, Parameters * param, lua_State * L) 
{
    if(luaL_loadstring(L, confstr) || lua_pcall(L, 0, 1, 0)) {
        /* This shall never happen unless the dump library is brokean */
        fastpm_raise(-1, "%s\n", lua_tostring(L, -1));
    }
    
    fastpm_log(-1, "Configuration %s\n", confstr);

    memset(param, 0, sizeof(*param));
    param->nc = read_integer(L, "nc");
    param->boxsize = read_number(L, "boxsize");

    param->time_step = read_array_number(L, "time_step", &param->n_time_step);

    param->aout = read_array_number(L, "output_redshifts", &param->n_aout);
    /* convert z to a */
    int i;
    for(i = 0; i < param->n_aout; i ++) {
        param->aout[i] = 1 / (param->aout[i] + 1);
    }
    param->random_seed = read_integer(L, "random_seed");

    param->omega_m = read_number(L, "omega_m");
    param->h = read_number(L, "h");
    param->sigma8 = read_number_opt(L, "sigma8", 0);

    param->pm_nc_factor = read_array_integer(L, "pm_nc_factor", &param->n_pm_nc_factor);
    param->change_pm = read_array_number(L, "change_pm", &param->n_change_pm);
    if(param->n_change_pm != param->n_pm_nc_factor) {
        fastpm_raise(1001, "Error: length of change_pm and pm_nc_factor mismatch.");
    }
    param->np_alloc_factor = read_number(L, "np_alloc_factor");

    // File Names and optional parameters realated
    param->read_grafic= read_string_opt(L, "read_grafic", NULL);
    param->read_noisek = read_string_opt(L, "read_noisek", NULL);
    param->read_noise = read_string_opt(L, "read_noise", NULL);
    param->read_runpbic = read_string_opt(L, "read_runpbic", NULL);

    /* compatible */
    param->read_runpbic = read_string_opt(L, "readic", param->read_runpbic);

    if(!param->read_runpbic) {
        param->read_powerspectrum = read_string_opt(L, "read_powerspectrum", NULL);
        param->read_powerspectrum = read_string_opt(L, "powerspectrum", param->read_powerspectrum);
    }

    param->write_powerspectrum = read_string_opt(L, "write_powerspectrum", NULL);
    param->write_powerspectrum = read_string_opt(L, "measure_power", param->write_powerspectrum);

    param->write_snapshot = read_string_opt(L, "write_snapshot", NULL);
    param->write_snapshot = read_string_opt(L, "snapshot", param->write_snapshot);
    param->write_snapshot_runpb = read_string_opt(L, "write_snapshot_runpb", NULL);
    param->write_noisek = read_string_opt(L, "write_noisek", NULL);
    param->write_noise = read_string_opt(L, "write_noise", NULL);

    struct enum_entry table[] = {
        {"cola", FORCE_MODE_COLA},
        {"pm", FORCE_MODE_PM},
        {NULL, -1},
    };

    param->force_mode = read_enum(L, "force_mode", table);
    
    if(param->force_mode == FORCE_MODE_PM) {
        param->cola_stdda = 1;
        param->enforce_broadband = 1;
    } else {
        param->cola_stdda = 0;
        param->enforce_broadband = 0;
    }
    param->cola_stdda = read_boolean_opt(L, "cola_stdda", param->cola_stdda);
    param->enforce_broadband = read_boolean_opt(L, "enforce_broadband", param->enforce_broadband);
    param->enforce_broadband_kmax = read_integer_opt(L, "enforce_broadband_kmax", 4);

    param->use_dx1_only = read_boolean_opt(L, "za", 0);
}

