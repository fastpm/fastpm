    #include <lua.h>
    #include <lauxlib.h>
    #include <lualib.h>
    #include <string.h>
    #include <stdlib.h>

    typedef struct LuaConfig LuaConfig;

    typedef struct LuaConfigRef LuaConfigRef;

    struct LuaConfigRef {
        char * name;
        void * data;
        struct LuaConfigRef * next;
    };

    struct LuaConfig {
        lua_State * L;
        char * error;
        LuaConfigRef head;
    };

    /* Helper function s */
    static int luaL_eval(lua_State * L, const char * string)
    {
        /* Evaluate string based on the global namespace at stack top */
        char * s = malloc(strlen(string) + 20);
        sprintf(s, "return %s", string);
        luaL_loadstring(L, s);
        free(s);
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

    LuaConfig * lua_config_new(const char * luastring)
    {
        LuaConfig * lc = malloc(sizeof(lc[0]));
        lua_State * L = luaL_newstate();
        luaL_openlibs(L);

        lc->L = L;

        lc->error = NULL;

        lc->head.next = NULL;
        if(luaL_eval(L, luastring)) {
            lc->error = _strdup(lua_tostring(L, -1));
        }
        return lc;
    }

    char *
    lua_config_parse(char * entrypoint, char * filename, int argc, char ** argv, char ** error)
    {

        char * confstr;

        extern int lua_open_runtime(lua_State * L);

        lua_State *L = luaL_newstate();
        luaL_openlibs(L);

        if(lua_open_runtime(L)) {
            *error = _strdup(lua_tostring(L, -1));
            goto fail;
        }

        lua_getglobal(L, entrypoint);
        char * real = filename; //realpath(filename, NULL);
        lua_pushstring(L, real);

        int i;

        for(i = 0; i < argc; i ++) {
            lua_pushstring(L, argv[i]);
        }

        if(lua_pcall(L, 1 + argc, 1, 0)) {
            *error = _strdup(lua_tostring(L, -1));
            goto fail;
        }
        confstr = _strdup(lua_tostring(L, -1));
        lua_pop(L, 1);

        lua_close(L);
        return confstr;

    fail:
        lua_close(L);
        return NULL;
    }
    static void
    fill_shape_and_strides(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
    {
        int k = 0;
        int istable = 0;
        while(1) {
            istable = lua_istable(lc->L, -1);
            if(istable) {
                shape[k] = luaL_len(lc->L, -1);
                k ++;
                lua_pushinteger(lc->L, 1);
                lua_gettable(lc->L, -2);
            } else {
                break;
            }
        }

        *ndim = k;
        *size = 1;
        for(k = 0; k < *ndim; k ++) {
            *size *= shape[k];
        }
        if (*ndim >= 1) {
            strides[*ndim - 1] = 1;
            for(k = *ndim - 2; k >= 0; k --) {
                strides[k] = strides[k + 1] * shape[k + 1];
            }
        }
        lua_pop(lc->L, *ndim);
    }

    void lua_config_free(LuaConfig * lc)
    {
        if(lc->error) free(lc->error);
        LuaConfigRef * entry = lc->head.next;
        while(entry) {
            LuaConfigRef * q = entry->next;
            free(entry->name);
            /* unref the pointer */
            free(entry->data);
            free(entry);
            entry = q;
        }
        lua_close(lc->L);
        free(lc);
    }
    const char * lua_config_error(LuaConfig * lc)
    {
        return lc->error;
    }
    static void *
    lua_config_cache_set(LuaConfig * lc, const char * name, void * data)
    {
        LuaConfigRef * entry = lc->head.next;
        while(entry) {
            if(0 == strcmp(entry->name, name)) {
                if(data != entry->data) {
                    /* unref the old pointer */
                    free(entry->data);
                    entry->data = data;
                }
                return data;
            }
            entry = entry->next;
        }
        if(entry == NULL) {
            entry = malloc(sizeof(entry[0]));
            entry->name = _strdup(name);
            entry->data = data;
            entry->next = lc->head.next;
            lc->head.next = entry;
        }
        return data;
    }
    static void *
    lua_config_cache_get(LuaConfig * lc, const char * name) {
        LuaConfigRef * entry = lc->head.next;
        while(entry) {
            if(0 == strcmp(entry->name, name)) {
                return entry->data;
            }
            entry = entry->next;
        }
        return NULL;
    }

#include <fastpm/libfastpm.h>
            int lua_config_has_nc(LuaConfig * lc)
            {
                luaL_eval(lc->L, "nc");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_nc(LuaConfig * lc)
            {
                luaL_eval(lc->L, "nc");
                double val = lua_tointeger(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_boxsize(LuaConfig * lc)
            {
                luaL_eval(lc->L, "boxsize");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_boxsize(LuaConfig * lc)
            {
                luaL_eval(lc->L, "boxsize");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_time_step(LuaConfig * lc)
            {
                luaL_eval(lc->L, "time_step");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_time_step(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_time_step(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_time_step(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_time_step_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "time_step");
            if(cached) return cached;

            luaL_eval(lc->L, "time_step");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_time_step(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "time_step_ndim", ndim);
            lua_config_cache_set(lc, "time_step_shape", shape);
            lua_config_cache_set(lc, "time_step_strides", strides);
            lua_config_cache_set(lc, "time_step_size", size);
            return array;
        }
        static void lua_config_cache_time_step(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_time_step_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "time_step", array);
            lua_config_cache_set(lc, "time_step_ndim", ndim);
            lua_config_cache_set(lc, "time_step_shape", shape);
            lua_config_cache_set(lc, "time_step_strides", strides);
            lua_config_cache_set(lc, "time_step_size", size);
        }
        double * lua_config_get_time_step(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "time_step")) {
                lua_config_cache_time_step(lc);
            }
            return (double* )lua_config_cache_get(lc, "time_step");
        }

        int lua_config_get_n_time_step(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "time_step")) {
                lua_config_cache_time_step(lc);
            }
            return *(int*)lua_config_cache_get(lc, "time_step_shape");
        }
        int lua_config_get_ndim_time_step(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "time_step")) {
                lua_config_cache_time_step(lc);
            }
            return *(int*)lua_config_cache_get(lc, "time_step_ndim");
        }
        int * lua_config_get_shape_time_step(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "time_step")) {
                lua_config_cache_time_step(lc);
            }
            return (int*)lua_config_cache_get(lc, "time_step_shape");
        }
            int lua_config_has_output_redshifts(LuaConfig * lc)
            {
                luaL_eval(lc->L, "output_redshifts");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_output_redshifts(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_output_redshifts(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_output_redshifts(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_output_redshifts_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "output_redshifts");
            if(cached) return cached;

            luaL_eval(lc->L, "output_redshifts");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_output_redshifts(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "output_redshifts_ndim", ndim);
            lua_config_cache_set(lc, "output_redshifts_shape", shape);
            lua_config_cache_set(lc, "output_redshifts_strides", strides);
            lua_config_cache_set(lc, "output_redshifts_size", size);
            return array;
        }
        static void lua_config_cache_output_redshifts(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_output_redshifts_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "output_redshifts", array);
            lua_config_cache_set(lc, "output_redshifts_ndim", ndim);
            lua_config_cache_set(lc, "output_redshifts_shape", shape);
            lua_config_cache_set(lc, "output_redshifts_strides", strides);
            lua_config_cache_set(lc, "output_redshifts_size", size);
        }
        double * lua_config_get_output_redshifts(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "output_redshifts")) {
                lua_config_cache_output_redshifts(lc);
            }
            return (double* )lua_config_cache_get(lc, "output_redshifts");
        }

        int lua_config_get_n_output_redshifts(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "output_redshifts")) {
                lua_config_cache_output_redshifts(lc);
            }
            return *(int*)lua_config_cache_get(lc, "output_redshifts_shape");
        }
        int lua_config_get_ndim_output_redshifts(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "output_redshifts")) {
                lua_config_cache_output_redshifts(lc);
            }
            return *(int*)lua_config_cache_get(lc, "output_redshifts_ndim");
        }
        int * lua_config_get_shape_output_redshifts(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "output_redshifts")) {
                lua_config_cache_output_redshifts(lc);
            }
            return (int*)lua_config_cache_get(lc, "output_redshifts_shape");
        }
            int lua_config_has_aout(LuaConfig * lc)
            {
                luaL_eval(lc->L, "aout");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_aout(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_aout(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_aout(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_aout_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "aout");
            if(cached) return cached;

            luaL_eval(lc->L, "aout");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_aout(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "aout_ndim", ndim);
            lua_config_cache_set(lc, "aout_shape", shape);
            lua_config_cache_set(lc, "aout_strides", strides);
            lua_config_cache_set(lc, "aout_size", size);
            return array;
        }
        static void lua_config_cache_aout(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_aout_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "aout", array);
            lua_config_cache_set(lc, "aout_ndim", ndim);
            lua_config_cache_set(lc, "aout_shape", shape);
            lua_config_cache_set(lc, "aout_strides", strides);
            lua_config_cache_set(lc, "aout_size", size);
        }
        double * lua_config_get_aout(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "aout")) {
                lua_config_cache_aout(lc);
            }
            return (double* )lua_config_cache_get(lc, "aout");
        }

        int lua_config_get_n_aout(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "aout")) {
                lua_config_cache_aout(lc);
            }
            return *(int*)lua_config_cache_get(lc, "aout_shape");
        }
        int lua_config_get_ndim_aout(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "aout")) {
                lua_config_cache_aout(lc);
            }
            return *(int*)lua_config_cache_get(lc, "aout_ndim");
        }
        int * lua_config_get_shape_aout(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "aout")) {
                lua_config_cache_aout(lc);
            }
            return (int*)lua_config_cache_get(lc, "aout_shape");
        }
            int lua_config_has_omega_m(LuaConfig * lc)
            {
                luaL_eval(lc->L, "omega_m");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_omega_m(LuaConfig * lc)
            {
                luaL_eval(lc->L, "omega_m");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_h(LuaConfig * lc)
            {
                luaL_eval(lc->L, "h");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_h(LuaConfig * lc)
            {
                luaL_eval(lc->L, "h");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pm_nc_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pm_nc_factor");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_pm_nc_factor(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_pm_nc_factor(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_pm_nc_factor(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_pm_nc_factor_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "pm_nc_factor");
            if(cached) return cached;

            luaL_eval(lc->L, "pm_nc_factor");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_pm_nc_factor(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "pm_nc_factor_ndim", ndim);
            lua_config_cache_set(lc, "pm_nc_factor_shape", shape);
            lua_config_cache_set(lc, "pm_nc_factor_strides", strides);
            lua_config_cache_set(lc, "pm_nc_factor_size", size);
            return array;
        }
        static void lua_config_cache_pm_nc_factor(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_pm_nc_factor_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "pm_nc_factor", array);
            lua_config_cache_set(lc, "pm_nc_factor_ndim", ndim);
            lua_config_cache_set(lc, "pm_nc_factor_shape", shape);
            lua_config_cache_set(lc, "pm_nc_factor_strides", strides);
            lua_config_cache_set(lc, "pm_nc_factor_size", size);
        }
        double * lua_config_get_pm_nc_factor(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "pm_nc_factor")) {
                lua_config_cache_pm_nc_factor(lc);
            }
            return (double* )lua_config_cache_get(lc, "pm_nc_factor");
        }

        int lua_config_get_n_pm_nc_factor(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "pm_nc_factor")) {
                lua_config_cache_pm_nc_factor(lc);
            }
            return *(int*)lua_config_cache_get(lc, "pm_nc_factor_shape");
        }
        int lua_config_get_ndim_pm_nc_factor(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "pm_nc_factor")) {
                lua_config_cache_pm_nc_factor(lc);
            }
            return *(int*)lua_config_cache_get(lc, "pm_nc_factor_ndim");
        }
        int * lua_config_get_shape_pm_nc_factor(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "pm_nc_factor")) {
                lua_config_cache_pm_nc_factor(lc);
            }
            return (int*)lua_config_cache_get(lc, "pm_nc_factor_shape");
        }
            int lua_config_has_np_alloc_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "np_alloc_factor");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_np_alloc_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "np_alloc_factor");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_compute_potential(LuaConfig * lc)
            {
                luaL_eval(lc->L, "compute_potential");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_compute_potential(LuaConfig * lc)
            {
                luaL_eval(lc->L, "compute_potential");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_m_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "m_ncdm");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_m_ncdm(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_m_ncdm(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_m_ncdm(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_m_ncdm_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "m_ncdm");
            if(cached) return cached;

            luaL_eval(lc->L, "m_ncdm");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_m_ncdm(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "m_ncdm_ndim", ndim);
            lua_config_cache_set(lc, "m_ncdm_shape", shape);
            lua_config_cache_set(lc, "m_ncdm_strides", strides);
            lua_config_cache_set(lc, "m_ncdm_size", size);
            return array;
        }
        static void lua_config_cache_m_ncdm(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_m_ncdm_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "m_ncdm", array);
            lua_config_cache_set(lc, "m_ncdm_ndim", ndim);
            lua_config_cache_set(lc, "m_ncdm_shape", shape);
            lua_config_cache_set(lc, "m_ncdm_strides", strides);
            lua_config_cache_set(lc, "m_ncdm_size", size);
        }
        double * lua_config_get_m_ncdm(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "m_ncdm")) {
                lua_config_cache_m_ncdm(lc);
            }
            return (double* )lua_config_cache_get(lc, "m_ncdm");
        }

        int lua_config_get_n_m_ncdm(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "m_ncdm")) {
                lua_config_cache_m_ncdm(lc);
            }
            return *(int*)lua_config_cache_get(lc, "m_ncdm_shape");
        }
        int lua_config_get_ndim_m_ncdm(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "m_ncdm")) {
                lua_config_cache_m_ncdm(lc);
            }
            return *(int*)lua_config_cache_get(lc, "m_ncdm_ndim");
        }
        int * lua_config_get_shape_m_ncdm(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "m_ncdm")) {
                lua_config_cache_m_ncdm(lc);
            }
            return (int*)lua_config_cache_get(lc, "m_ncdm_shape");
        }
            int lua_config_has_every_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "every_ncdm");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_every_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "every_ncdm");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_n_shell(LuaConfig * lc)
            {
                luaL_eval(lc->L, "n_shell");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_n_shell(LuaConfig * lc)
            {
                luaL_eval(lc->L, "n_shell");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_n_side(LuaConfig * lc)
            {
                luaL_eval(lc->L, "n_side");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_n_side(LuaConfig * lc)
            {
                luaL_eval(lc->L, "n_side");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lvk(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lvk");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_lvk(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lvk");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_ncdm_sphere_scheme(LuaConfig * lc)
            {
                luaL_eval(lc->L, "ncdm_sphere_scheme");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_ncdm_sphere_scheme(LuaConfig * lc)
            {
                luaL_eval(lc->L, "ncdm_sphere_scheme");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("fibonacci", str)) {
                    val = FASTPM_NCDM_SPHERE_FIBONACCI;
                }
                if (0 == strcmp("healpix", str)) {
                    val = FASTPM_NCDM_SPHERE_HEALPIX;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_painter_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "painter_type");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_painter_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "painter_type");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("cic", str)) {
                    val = FASTPM_PAINTER_CIC;
                }
                if (0 == strcmp("lanczos", str)) {
                    val = FASTPM_PAINTER_LANCZOS;
                }
                if (0 == strcmp("linear", str)) {
                    val = FASTPM_PAINTER_LINEAR;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_painter_support(LuaConfig * lc)
            {
                luaL_eval(lc->L, "painter_support");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_painter_support(LuaConfig * lc)
            {
                luaL_eval(lc->L, "painter_support");
                double val = lua_tointeger(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_force_mode(LuaConfig * lc)
            {
                luaL_eval(lc->L, "force_mode");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_force_mode(LuaConfig * lc)
            {
                luaL_eval(lc->L, "force_mode");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("cola", str)) {
                    val = FASTPM_FORCE_COLA;
                }
                if (0 == strcmp("zola", str)) {
                    val = FASTPM_FORCE_FASTPM;
                }
                if (0 == strcmp("pm", str)) {
                    val = FASTPM_FORCE_PM;
                }
                if (0 == strcmp("fastpm", str)) {
                    val = FASTPM_FORCE_FASTPM;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_enforce_broadband_kmax(LuaConfig * lc)
            {
                luaL_eval(lc->L, "enforce_broadband_kmax");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_enforce_broadband_kmax(LuaConfig * lc)
            {
                luaL_eval(lc->L, "enforce_broadband_kmax");
                double val = lua_tointeger(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_f_nl_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "f_nl_type");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_f_nl_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "f_nl_type");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("none", str)) {
                    val = FASTPM_FNL_NONE;
                }
                if (0 == strcmp("local", str)) {
                    val = FASTPM_FNL_LOCAL;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_f_nl(LuaConfig * lc)
            {
                luaL_eval(lc->L, "f_nl");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_f_nl(LuaConfig * lc)
            {
                luaL_eval(lc->L, "f_nl");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_kmax_primordial_over_knyquist(LuaConfig * lc)
            {
                luaL_eval(lc->L, "kmax_primordial_over_knyquist");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_kmax_primordial_over_knyquist(LuaConfig * lc)
            {
                luaL_eval(lc->L, "kmax_primordial_over_knyquist");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_scalar_amp(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_amp");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_scalar_amp(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_amp");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_scalar_pivot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_pivot");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_scalar_pivot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_pivot");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_scalar_spectral_index(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_spectral_index");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_scalar_spectral_index(LuaConfig * lc)
            {
                luaL_eval(lc->L, "scalar_spectral_index");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_read_lineark(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_lineark");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_lineark(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_lineark");
                if(cached) return cached;
                luaL_eval(lc->L, "read_lineark");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_lineark", _strdup(val));
                return NULL;
            }
            int lua_config_has_read_powerspectrum(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_powerspectrum");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_powerspectrum(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_powerspectrum");
                if(cached) return cached;
                luaL_eval(lc->L, "read_powerspectrum");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_powerspectrum", _strdup(val));
                return NULL;
            }
            int lua_config_has_linear_density_redshift(LuaConfig * lc)
            {
                luaL_eval(lc->L, "linear_density_redshift");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_linear_density_redshift(LuaConfig * lc)
            {
                luaL_eval(lc->L, "linear_density_redshift");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_read_lineark_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_lineark_ncdm");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_lineark_ncdm(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_lineark_ncdm");
                if(cached) return cached;
                luaL_eval(lc->L, "read_lineark_ncdm");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_lineark_ncdm", _strdup(val));
                return NULL;
            }
            int lua_config_has_read_powerspectrum_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_powerspectrum_ncdm");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_powerspectrum_ncdm(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_powerspectrum_ncdm");
                if(cached) return cached;
                luaL_eval(lc->L, "read_powerspectrum_ncdm");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_powerspectrum_ncdm", _strdup(val));
                return NULL;
            }
            int lua_config_has_linear_density_redshift_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "linear_density_redshift_ncdm");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_linear_density_redshift_ncdm(LuaConfig * lc)
            {
                luaL_eval(lc->L, "linear_density_redshift_ncdm");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_read_grafic(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_grafic");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_grafic(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_grafic");
                if(cached) return cached;
                luaL_eval(lc->L, "read_grafic");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_grafic", _strdup(val));
                return NULL;
            }
            int lua_config_has_read_runpbic(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_runpbic");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_runpbic(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_runpbic");
                if(cached) return cached;
                luaL_eval(lc->L, "read_runpbic");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_runpbic", _strdup(val));
                return NULL;
            }
            int lua_config_has_read_whitenoisek(LuaConfig * lc)
            {
                luaL_eval(lc->L, "read_whitenoisek");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_read_whitenoisek(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "read_whitenoisek");
                if(cached) return cached;
                luaL_eval(lc->L, "read_whitenoisek");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "read_whitenoisek", _strdup(val));
                return NULL;
            }
            int lua_config_has_sigma8(LuaConfig * lc)
            {
                luaL_eval(lc->L, "sigma8");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_sigma8(LuaConfig * lc)
            {
                luaL_eval(lc->L, "sigma8");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_random_seed(LuaConfig * lc)
            {
                luaL_eval(lc->L, "random_seed");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_random_seed(LuaConfig * lc)
            {
                luaL_eval(lc->L, "random_seed");
                double val = lua_tointeger(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_shift(LuaConfig * lc)
            {
                luaL_eval(lc->L, "shift");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_shift(LuaConfig * lc)
            {
                luaL_eval(lc->L, "shift");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_inverted_ic(LuaConfig * lc)
            {
                luaL_eval(lc->L, "inverted_ic");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_inverted_ic(LuaConfig * lc)
            {
                luaL_eval(lc->L, "inverted_ic");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_remove_cosmic_variance(LuaConfig * lc)
            {
                luaL_eval(lc->L, "remove_cosmic_variance");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_remove_cosmic_variance(LuaConfig * lc)
            {
                luaL_eval(lc->L, "remove_cosmic_variance");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_write_lineark(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_lineark");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_lineark(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_lineark");
                if(cached) return cached;
                luaL_eval(lc->L, "write_lineark");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_lineark", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_whitenoisek(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_whitenoisek");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_whitenoisek(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_whitenoisek");
                if(cached) return cached;
                luaL_eval(lc->L, "write_whitenoisek");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_whitenoisek", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_runpbic(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_runpbic");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_runpbic(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_runpbic");
                if(cached) return cached;
                luaL_eval(lc->L, "write_runpbic");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_runpbic", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_powerspectrum(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_powerspectrum");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_powerspectrum(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_powerspectrum");
                if(cached) return cached;
                luaL_eval(lc->L, "write_powerspectrum");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_powerspectrum", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_snapshot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_snapshot");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_snapshot(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_snapshot");
                if(cached) return cached;
                luaL_eval(lc->L, "write_snapshot");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_snapshot", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_nonlineark(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_nonlineark");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_nonlineark(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_nonlineark");
                if(cached) return cached;
                luaL_eval(lc->L, "write_nonlineark");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_nonlineark", _strdup(val));
                return NULL;
            }
            int lua_config_has_write_runpb_snapshot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_runpb_snapshot");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_runpb_snapshot(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_runpb_snapshot");
                if(cached) return cached;
                luaL_eval(lc->L, "write_runpb_snapshot");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_runpb_snapshot", _strdup(val));
                return NULL;
            }
            int lua_config_has_particle_fraction(LuaConfig * lc)
            {
                luaL_eval(lc->L, "particle_fraction");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_particle_fraction(LuaConfig * lc)
            {
                luaL_eval(lc->L, "particle_fraction");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_sort_snapshot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "sort_snapshot");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_sort_snapshot(LuaConfig * lc)
            {
                luaL_eval(lc->L, "sort_snapshot");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_write_fof(LuaConfig * lc)
            {
                luaL_eval(lc->L, "write_fof");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_write_fof(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "write_fof");
                if(cached) return cached;
                luaL_eval(lc->L, "write_fof");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "write_fof", _strdup(val));
                return NULL;
            }
            int lua_config_has_fof_linkinglength(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_linkinglength");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_fof_linkinglength(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_linkinglength");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_fof_nmin(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_nmin");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_fof_nmin(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_nmin");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_fof_kdtree_thresh(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_kdtree_thresh");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_fof_kdtree_thresh(LuaConfig * lc)
            {
                luaL_eval(lc->L, "fof_kdtree_thresh");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_amin(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_amin");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_lc_amin(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_amin");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_amax(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_amax");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_lc_amax(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_amax");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_write_usmesh(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_write_usmesh");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_lc_write_usmesh(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "lc_write_usmesh");
                if(cached) return cached;
                luaL_eval(lc->L, "lc_write_usmesh");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "lc_write_usmesh", _strdup(val));
                return NULL;
            }
            int lua_config_has_lc_usmesh_alloc_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_usmesh_alloc_factor");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_lc_usmesh_alloc_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_usmesh_alloc_factor");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_usmesh_fof_padding(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_usmesh_fof_padding");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_lc_usmesh_fof_padding(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_usmesh_fof_padding");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_usmesh_tiles(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_usmesh_tiles");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_lc_usmesh_tiles(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_lc_usmesh_tiles(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_lc_usmesh_tiles(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_lc_usmesh_tiles_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "lc_usmesh_tiles");
            if(cached) return cached;

            luaL_eval(lc->L, "lc_usmesh_tiles");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_lc_usmesh_tiles(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "lc_usmesh_tiles_ndim", ndim);
            lua_config_cache_set(lc, "lc_usmesh_tiles_shape", shape);
            lua_config_cache_set(lc, "lc_usmesh_tiles_strides", strides);
            lua_config_cache_set(lc, "lc_usmesh_tiles_size", size);
            return array;
        }
        static void lua_config_cache_lc_usmesh_tiles(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_lc_usmesh_tiles_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "lc_usmesh_tiles", array);
            lua_config_cache_set(lc, "lc_usmesh_tiles_ndim", ndim);
            lua_config_cache_set(lc, "lc_usmesh_tiles_shape", shape);
            lua_config_cache_set(lc, "lc_usmesh_tiles_strides", strides);
            lua_config_cache_set(lc, "lc_usmesh_tiles_size", size);
        }
        double * lua_config_get_lc_usmesh_tiles(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_usmesh_tiles")) {
                lua_config_cache_lc_usmesh_tiles(lc);
            }
            return (double* )lua_config_cache_get(lc, "lc_usmesh_tiles");
        }

        int lua_config_get_n_lc_usmesh_tiles(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_usmesh_tiles")) {
                lua_config_cache_lc_usmesh_tiles(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_usmesh_tiles_shape");
        }
        int lua_config_get_ndim_lc_usmesh_tiles(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_usmesh_tiles")) {
                lua_config_cache_lc_usmesh_tiles(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_usmesh_tiles_ndim");
        }
        int * lua_config_get_shape_lc_usmesh_tiles(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_usmesh_tiles")) {
                lua_config_cache_lc_usmesh_tiles(lc);
            }
            return (int*)lua_config_cache_get(lc, "lc_usmesh_tiles_shape");
        }
            int lua_config_has_dh_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "dh_factor");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_dh_factor(LuaConfig * lc)
            {
                luaL_eval(lc->L, "dh_factor");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_fov(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_fov");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_lc_fov(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_fov");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_lc_octants(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_octants");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_lc_octants(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_lc_octants(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_lc_octants(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_lc_octants_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "lc_octants");
            if(cached) return cached;

            luaL_eval(lc->L, "lc_octants");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_lc_octants(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "lc_octants_ndim", ndim);
            lua_config_cache_set(lc, "lc_octants_shape", shape);
            lua_config_cache_set(lc, "lc_octants_strides", strides);
            lua_config_cache_set(lc, "lc_octants_size", size);
            return array;
        }
        static void lua_config_cache_lc_octants(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_lc_octants_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "lc_octants", array);
            lua_config_cache_set(lc, "lc_octants_ndim", ndim);
            lua_config_cache_set(lc, "lc_octants_shape", shape);
            lua_config_cache_set(lc, "lc_octants_strides", strides);
            lua_config_cache_set(lc, "lc_octants_size", size);
        }
        double * lua_config_get_lc_octants(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_octants")) {
                lua_config_cache_lc_octants(lc);
            }
            return (double* )lua_config_cache_get(lc, "lc_octants");
        }

        int lua_config_get_n_lc_octants(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_octants")) {
                lua_config_cache_lc_octants(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_octants_shape");
        }
        int lua_config_get_ndim_lc_octants(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_octants")) {
                lua_config_cache_lc_octants(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_octants_ndim");
        }
        int * lua_config_get_shape_lc_octants(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_octants")) {
                lua_config_cache_lc_octants(lc);
            }
            return (int*)lua_config_cache_get(lc, "lc_octants_shape");
        }
            int lua_config_has_lc_glmatrix(LuaConfig * lc)
            {
                luaL_eval(lc->L, "lc_glmatrix");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_lc_glmatrix(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_lc_glmatrix(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_lc_glmatrix(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_lc_glmatrix_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "lc_glmatrix");
            if(cached) return cached;

            luaL_eval(lc->L, "lc_glmatrix");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_lc_glmatrix(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "lc_glmatrix_ndim", ndim);
            lua_config_cache_set(lc, "lc_glmatrix_shape", shape);
            lua_config_cache_set(lc, "lc_glmatrix_strides", strides);
            lua_config_cache_set(lc, "lc_glmatrix_size", size);
            return array;
        }
        static void lua_config_cache_lc_glmatrix(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_lc_glmatrix_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "lc_glmatrix", array);
            lua_config_cache_set(lc, "lc_glmatrix_ndim", ndim);
            lua_config_cache_set(lc, "lc_glmatrix_shape", shape);
            lua_config_cache_set(lc, "lc_glmatrix_strides", strides);
            lua_config_cache_set(lc, "lc_glmatrix_size", size);
        }
        double * lua_config_get_lc_glmatrix(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_glmatrix")) {
                lua_config_cache_lc_glmatrix(lc);
            }
            return (double* )lua_config_cache_get(lc, "lc_glmatrix");
        }

        int lua_config_get_n_lc_glmatrix(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_glmatrix")) {
                lua_config_cache_lc_glmatrix(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_glmatrix_shape");
        }
        int lua_config_get_ndim_lc_glmatrix(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_glmatrix")) {
                lua_config_cache_lc_glmatrix(lc);
            }
            return *(int*)lua_config_cache_get(lc, "lc_glmatrix_ndim");
        }
        int * lua_config_get_shape_lc_glmatrix(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "lc_glmatrix")) {
                lua_config_cache_lc_glmatrix(lc);
            }
            return (int*)lua_config_cache_get(lc, "lc_glmatrix_shape");
        }
            int lua_config_has_za(LuaConfig * lc)
            {
                luaL_eval(lc->L, "za");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_za(LuaConfig * lc)
            {
                luaL_eval(lc->L, "za");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_kernel_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "kernel_type");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_kernel_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "kernel_type");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("3_4", str)) {
                    val = FASTPM_KERNEL_3_4;
                }
                if (0 == strcmp("3_2", str)) {
                    val = FASTPM_KERNEL_3_2;
                }
                if (0 == strcmp("5_4", str)) {
                    val = FASTPM_KERNEL_5_4;
                }
                if (0 == strcmp("eastwood", str)) {
                    val = FASTPM_KERNEL_EASTWOOD;
                }
                if (0 == strcmp("naive", str)) {
                    val = FASTPM_KERNEL_NAIVE;
                }
                if (0 == strcmp("gadget", str)) {
                    val = FASTPM_KERNEL_GADGET;
                }
                if (0 == strcmp("1_4", str)) {
                    val = FASTPM_KERNEL_1_4;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_force_softening_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "force_softening_type");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_force_softening_type(LuaConfig * lc)
            {
                luaL_eval(lc->L, "force_softening_type");
                const char * str = lua_tostring(lc->L, -1);
                int val;
                if (0 == strcmp("gadget_long_range", str)) {
                    val = FASTPM_SOFTENING_GADGET_LONG_RANGE;
                }
                if (0 == strcmp("twothird", str)) {
                    val = FASTPM_SOFTENING_TWO_THIRD;
                }
                if (0 == strcmp("gaussian36", str)) {
                    val = FASTPM_SOFTENING_GAUSSIAN36;
                }
                if (0 == strcmp("gaussian", str)) {
                    val = FASTPM_SOFTENING_GAUSSIAN;
                }
                if (0 == strcmp("none", str)) {
                    val = FASTPM_SOFTENING_NONE;
                }
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_constraints(LuaConfig * lc)
            {
                luaL_eval(lc->L, "constraints");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_constraints(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_constraints(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_constraints(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_constraints_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "constraints");
            if(cached) return cached;

            luaL_eval(lc->L, "constraints");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_constraints(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "constraints_ndim", ndim);
            lua_config_cache_set(lc, "constraints_shape", shape);
            lua_config_cache_set(lc, "constraints_strides", strides);
            lua_config_cache_set(lc, "constraints_size", size);
            return array;
        }
        static void lua_config_cache_constraints(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_constraints_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "constraints", array);
            lua_config_cache_set(lc, "constraints_ndim", ndim);
            lua_config_cache_set(lc, "constraints_shape", shape);
            lua_config_cache_set(lc, "constraints_strides", strides);
            lua_config_cache_set(lc, "constraints_size", size);
        }
        double * lua_config_get_constraints(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "constraints")) {
                lua_config_cache_constraints(lc);
            }
            return (double* )lua_config_cache_get(lc, "constraints");
        }

        int lua_config_get_n_constraints(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "constraints")) {
                lua_config_cache_constraints(lc);
            }
            return *(int*)lua_config_cache_get(lc, "constraints_shape");
        }
        int lua_config_get_ndim_constraints(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "constraints")) {
                lua_config_cache_constraints(lc);
            }
            return *(int*)lua_config_cache_get(lc, "constraints_ndim");
        }
        int * lua_config_get_shape_constraints(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "constraints")) {
                lua_config_cache_constraints(lc);
            }
            return (int*)lua_config_cache_get(lc, "constraints_shape");
        }
            int lua_config_has_set_mode_method(LuaConfig * lc)
            {
                luaL_eval(lc->L, "set_mode_method");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            const char * lua_config_get_set_mode_method(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "set_mode_method");
                if(cached) return cached;
                luaL_eval(lc->L, "set_mode_method");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "set_mode_method", _strdup(val));
                return NULL;
            }
            int lua_config_has_set_mode(LuaConfig * lc)
            {
                luaL_eval(lc->L, "set_mode");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }

        static int
        fill_data_set_mode(LuaConfig * lc, double * array, int ndim, int * shape, int * strides);

        static int
        fill_data_set_mode(LuaConfig * lc, double * array, int ndim, int * shape, int * strides)
        {
            if(ndim < 0) {
                return 1;
            }
            if(lua_istable(lc->L, -1)) {
                int i;
                int len = luaL_len(lc->L, -1);
                if(shape[0] != len) return 1;
                for(i = 0; i < len; i ++) {
                    lua_pushinteger(lc->L, i + 1);
                    lua_gettable(lc->L, -2);
                    if(0 != fill_data_set_mode(lc,
                            array + strides[0] * i,
                            ndim - 1,
                            shape + 1,
                            strides + 1)) {
                        lua_pop(lc->L, 1);
                        return 1;
                    }
                    lua_pop(lc->L, 1);
                }
            } else {
                double x = lua_tonumber(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        double * lua_config_get_set_mode_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            double * cached = (double *) lua_config_cache_get(lc, "set_mode");
            if(cached) return cached;

            luaL_eval(lc->L, "set_mode");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            double * array = (double *) malloc(sizeof(double) * *size);
            if(0 != fill_data_set_mode(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "set_mode_ndim", ndim);
            lua_config_cache_set(lc, "set_mode_shape", shape);
            lua_config_cache_set(lc, "set_mode_strides", strides);
            lua_config_cache_set(lc, "set_mode_size", size);
            return array;
        }
        static void lua_config_cache_set_mode(LuaConfig * lc)
        {
            double * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = lua_config_get_set_mode_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "set_mode", array);
            lua_config_cache_set(lc, "set_mode_ndim", ndim);
            lua_config_cache_set(lc, "set_mode_shape", shape);
            lua_config_cache_set(lc, "set_mode_strides", strides);
            lua_config_cache_set(lc, "set_mode_size", size);
        }
        double * lua_config_get_set_mode(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "set_mode")) {
                lua_config_cache_set_mode(lc);
            }
            return (double* )lua_config_cache_get(lc, "set_mode");
        }

        int lua_config_get_n_set_mode(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "set_mode")) {
                lua_config_cache_set_mode(lc);
            }
            return *(int*)lua_config_cache_get(lc, "set_mode_shape");
        }
        int lua_config_get_ndim_set_mode(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "set_mode")) {
                lua_config_cache_set_mode(lc);
            }
            return *(int*)lua_config_cache_get(lc, "set_mode_ndim");
        }
        int * lua_config_get_shape_set_mode(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "set_mode")) {
                lua_config_cache_set_mode(lc);
            }
            return (int*)lua_config_cache_get(lc, "set_mode_shape");
        }
            int lua_config_has_pgdc(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            int lua_config_get_pgdc(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pgdc_alpha0(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_alpha0");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_pgdc_alpha0(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_alpha0");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pgdc_A(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_A");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_pgdc_A(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_A");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pgdc_B(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_B");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_pgdc_B(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_B");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pgdc_kl(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_kl");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_pgdc_kl(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_kl");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
            int lua_config_has_pgdc_ks(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_ks");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
            double lua_config_get_pgdc_ks(LuaConfig * lc)
            {
                luaL_eval(lc->L, "pgdc_ks");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
