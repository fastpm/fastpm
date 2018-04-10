local dump = require("lua-runtime-dump")

local _NAME = ... or 'main'

config = {}

local function Schema()
    local self = {}

    local names = {}

    function self.declare(options)
        local name = options.name
        local type = options.type
        local required = options.required or false
        local default = options.default
        local choices = options.choices
        local action = options.action
        local help = options.help

        local subtype = options.type:match('array:(%a+)')
        if subtype ~= nil then
            type = 'array'
        end
        self[name] = {type=type, subtype=subtype, required=required, default=default, choices=choice, action=action, help=help}
        names[#names + 1] = name
    end

    function self.print()
        for i,k in pairs(names) do
            v = self[k]
            local header = string.format('name: %s ', k)
            local body = ''
            for kk,vv in pairs(v) do
                body = string.format('%s %s:%s', body, kk,vv)
            end
            print(string.format('%s - %s', header, body))
        end
    end

    function self.format_help()
        local template = [[ %-32s : %10s (%10s) : %s %s]] .. '\n'
        local result = ''
        for i, name in pairs(names) do
            local choices = ''
            if self[name].choices ~= nil then
                choices = ' Valid choices are : '
                for ch, v in pairs(self[name].choices) do
                    choices = choices .. ch .. ', '
                end
            end
            local line = string.format(template, name, self[name].type, tostring(self[name].default), self[name].help, choices)
            result = result .. line
        end
        return result
    end

    local function exists(filename)
        local file = io.open(filename, 'r')
        if file == nil then
            return true
            --return false
        end
        file:close()
        return true
    end

    local function check_array(name, entry, value)
        if type(value) == 'table' then
            for i, v in pairs(value) do
                check_array(name, entry, v)
            end
        else
            local subtype = entry.subtype
            if subtype == 'int' then
                subtype = 'number'
            end
            if type(value) ~= subtype then
                error(string.format("entry %d of key `%s` is of type `%s`, but `%s` is expected",
                    i, name, type(value), subtype))
            end
        end
    end

    local function check_int(name, entry, value)
        if type(value) ~= 'number' then
            error(string.format("key `%s` is of type `%s`, but `%s` is expected",
                name, type(value), entry.type))
        end
    end

    local function check_file(name, entry, value)
        if not exists(value) then
            error(string.format("file `%s' is not readable.", value))
        end
    end

    local function check_enum(name, entry, value)
        if entry.choices[value] == nil then
            local s = ''
            for v,_ in pairs(entry.choices) do
                s = s .. string.format('`%s`', v) .. ' '
            end
            error(string.format("value `%s` of key `%s` is not one of %s",
                value, name, s))
        end
    end

    local function check_type(name, entry, value)
        if entry.type == 'file' then
            -- File?
            check_file(name, entry, value)
        elseif entry.type == 'array' then
            -- Array?
            check_array(name, entry, value)
        elseif entry.type == 'enum' then
            -- Enum?
            check_enum(name, entry, value)
        elseif entry.type == 'int' then
            check_int(name, entry, value)
        elseif type(value) ~= entry.type then
            -- Anything else
            error(string.format("key `%s` is of type `%s`, but `%s` is expected",
                name, type(value), entry.type))
        end
    end

    local function bind_one(name, entry, value)
        if value == nil then
            if entry.required then
                error(string.format("`%s` is required but undefined.", name))
            else
                value = entry.default
            end
        end
        -- default could be nil
        -- nil matches any type requirements
        if value ~= nil then
            check_type(name, entry, value)
        end

        if entry.action ~= nil then
            entry.action(value)
        end
        return value
    end

    function self.bind(namespace)
        result = {}
        for i,k in pairs(names) do
            result[k] = bind_one(k, self[k], namespace[k])
        end
        return result
    end

    function self.pairs()
        local i
        local function next(self, name)
            if name == nil then
                i = 1
            end
            name = names[i]
            if name == nil then
                return nil
            end
            i = i + 1
            return name, self[name]
        end
        return next, self, nil
    end
    return self
end

config.Schema = Schema

local function rtrim(s)
    s = string.gsub(s, "^(.-)[ ]*$", "%1")
    return s
end

local function visit_string(name, entry, mode)
    if mode == 'h' then
        return [[
            const char * @PREFIX@_get_@name@(LuaConfig * lc);
        ]]
    else
        return [[
            const char * @PREFIX@_get_@name@(LuaConfig * lc)
            {
                char * cached = (char*) lua_config_cache_get(lc, "@name@");
                if(cached) return cached;
                luaL_eval(lc->L, "@name@");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return lua_config_cache_set(lc, "@name@", _strdup(val));
                return NULL;
            }
        ]]
    end
end

local function visit_number(name, entry, mode)
    if mode == 'h' then
        return [[
            double @PREFIX@_get_@name@(LuaConfig * lc);
        ]]
    else
        return [[
            double @PREFIX@_get_@name@(LuaConfig * lc)
            {
                luaL_eval(lc->L, "@name@");
                double val = lua_tonumber(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
        ]]
    end
end

local function visit_int(name, entry, mode)
    if mode == 'h' then
        return [[
            int @PREFIX@_get_@name@(LuaConfig * lc);
        ]]
    else
        return [[
            int @PREFIX@_get_@name@(LuaConfig * lc)
            {
                luaL_eval(lc->L, "@name@");
                double val = lua_tointeger(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
        ]]
    end
end

local function visit_array(name, entry, mode)
    local subctypetable = {
        number = 'double',
        int = 'int',
    }
    local convertertable = {
        number = 'lua_tonumber',
        int = 'lua_tointeger',
    }
    local function process(s)
        s = string.gsub(s, '@CTYPE@', subctypetable[entry.subtype])
        s = string.gsub(s, '@CONVERTOR@', convertertable[entry.subtype])
        return s
    end

    if mode == 'h' then
        return process( [[
            @CTYPE@ * @PREFIX@_get_@name@_full(LuaConfig * lc, int * size);
            @CTYPE@ * @PREFIX@_get_@name@(LuaConfig * lc);
            int @PREFIX@_get_n_@name@(LuaConfig * lc);
            int @PREFIX@_get_ndim_@name@(LuaConfig * lc);
            int * @PREFIX@_get_shape_@name@(LuaConfig * lc);
        ]])
    else
        return process([[

        static int
        fill_data_@name@(LuaConfig * lc, @CTYPE@ * array, int ndim, int * shape, int * strides);

        static int
        fill_data_@name@(LuaConfig * lc, @CTYPE@ * array, int ndim, int * shape, int * strides)
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
                    if(0 != fill_data_@name@(lc,
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
                @CTYPE@ x = @CONVERTOR@(lc->L, -1);
                array[0] = x;
            }
            return 0;
        }


        @CTYPE@ * @PREFIX@_get_@name@_full(LuaConfig * lc, int * ndim, int * shape, int * strides, int * size)
        {
            @CTYPE@ * cached = (@CTYPE@ *) lua_config_cache_get(lc, "@name@");
            if(cached) return cached;

            luaL_eval(lc->L, "@name@");
            if(lua_isnil(lc->L, -1)) {
                lua_pop(lc->L, 1);
                *size = 0;
                return NULL;
            };
            fill_shape_and_strides(lc, ndim, shape, strides, size);
            @CTYPE@ * array = (@CTYPE@ *) malloc(sizeof(@CTYPE@) * *size);
            if(0 != fill_data_@name@(lc, array, *ndim, shape, strides)) {
                free(array);
                lua_pop(lc->L,1);
                return NULL;
            }
            lua_pop(lc->L,1);
            lua_config_cache_set(lc, "@name@_ndim", ndim);
            lua_config_cache_set(lc, "@name@_shape", shape);
            lua_config_cache_set(lc, "@name@_strides", strides);
            lua_config_cache_set(lc, "@name@_size", size);
            return array;
        }
        static void @PREFIX@_cache_@name@(LuaConfig * lc)
        {
            @CTYPE@ * array = NULL;
            int * ndim = (int*) malloc(sizeof(int));
            int * shape = (int*) malloc(sizeof(int) * 32);
            int * strides = (int*) malloc(sizeof(int) * 32);
            int * size = (int*) malloc(sizeof(int));

            array = @PREFIX@_get_@name@_full(lc, ndim, shape, strides, size);

            lua_config_cache_set(lc, "@name@", array);
            lua_config_cache_set(lc, "@name@_ndim", ndim);
            lua_config_cache_set(lc, "@name@_shape", shape);
            lua_config_cache_set(lc, "@name@_strides", strides);
            lua_config_cache_set(lc, "@name@_size", size);
        }
        @CTYPE@ * @PREFIX@_get_@name@(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "@name@")) {
                @PREFIX@_cache_@name@(lc);
            }
            return (@CTYPE@* )lua_config_cache_get(lc, "@name@");
        }

        int @PREFIX@_get_n_@name@(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "@name@")) {
                @PREFIX@_cache_@name@(lc);
            }
            return *(int*)lua_config_cache_get(lc, "@name@_shape");
        }
        int @PREFIX@_get_ndim_@name@(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "@name@")) {
                @PREFIX@_cache_@name@(lc);
            }
            return *(int*)lua_config_cache_get(lc, "@name@_ndim");
        }
        int * @PREFIX@_get_shape_@name@(LuaConfig * lc)
        {
            if(!lua_config_cache_get(lc, "@name@")) {
                @PREFIX@_cache_@name@(lc);
            }
            return (int*)lua_config_cache_get(lc, "@name@_shape");
        }
        ]])
    end
end

local function visit_enum(name, entry, mode)
    if mode == 'h' then
        return [[
            int @PREFIX@_get_@name@(LuaConfig * lc);
        ]]
    else
        def = ''
        for k, v in pairs(entry.choices) do
            local e = [[
                if (0 == strcmp("%s", str)) {
                    val = %s;
                }
            ]]
            e = rtrim(string.format(e, k, v))
            def = def .. e
        end
        local a = [[
            int @PREFIX@_get_@name@(LuaConfig * lc)
            {
                luaL_eval(lc->L, "@name@");
                const char * str = lua_tostring(lc->L, -1);
                int val;
            ]]
        local b = [[
                lua_pop(lc->L, 1);
                return val;
            }
        ]]
        return rtrim(a) .. rtrim(def) .. rtrim(b)
    end
end

local function visit_boolean(name, entry, mode)
    if mode == 'h' then
        return [[
            int @PREFIX@_get_@name@(LuaConfig * lc);
        ]]
    else
        return [[
            int @PREFIX@_get_@name@(LuaConfig * lc)
            {
                luaL_eval(lc->L, "@name@");
                double val = lua_toboolean(lc->L, -1);
                lua_pop(lc->L, 1);
                return val;
            }
        ]]
    end
end
local function visit_any(name, entry, mode)
    if mode == 'h' then
        return [[
            int @PREFIX@_has_@name@(LuaConfig * lc);
        ]]
    else
        return [[
            int @PREFIX@_has_@name@(LuaConfig * lc)
            {
                luaL_eval(lc->L, "@name@");
                int isnil = lua_isnil(lc->L, -1);
                lua_pop(lc->L, 1);
                return !isnil;
            }
        ]]
    end
end

local function visit(name, entry, mode)
    if entry.type == 'string' then
        return visit_string(name, entry, mode)
    elseif entry.type == 'int' then
        return visit_int(name, entry, mode)
    elseif entry.type == 'boolean' then
        return visit_boolean(name, entry, mode)
    elseif entry.type == 'enum' then
        return visit_enum(name, entry, mode)
    elseif entry.type == 'array' then
        return visit_array(name, entry, mode)
    elseif entry.type == 'file' then
        return visit_string(name, entry, mode)
    elseif entry.type == 'number' then
        return visit_number(name, entry, mode)
    else
        error(name)
    end
end

function config.compile(schema, opt)
    local preample_h = [[
        typedef struct LuaConfig LuaConfig;

        void lua_config_free(LuaConfig * lc);
        LuaConfig * lua_config_new(const char * luastring);
        const char * lua_config_error(LuaConfig * lc);
        /* call global lua function with filename and arg as arguments */
        char * lua_config_parse(char * entrypoint, char * filename, int argc, char ** argv, char ** error);
    ]]

    local preample_c = [[
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

    ]]

    local stream_h = rtrim(preample_h)
    local stream_c = rtrim(preample_c)

    local prefix = opt.prefix or 'lua_config'

    for i, fn in pairs(opt.global_headers or {}) do
        local include = "#include <" .. fn .. ">\n"
        stream_c = stream_c .. include
    end
    for i, fn in pairs(opt.local_headers or {}) do
        local include = "#include \"" .. fn .. "\"\n"
        stream_c = stream_c .. include
    end

    for name, entry in schema.pairs() do
        local h1 = rtrim(visit_any(name, entry, 'h'))
        local c1 = rtrim(visit_any(name, entry, 'c'))
        local h = rtrim(visit(name, entry, 'h'))
        local c = rtrim(visit(name, entry, 'c'))

        h = string.gsub(h1 .. h, '@PREFIX@', prefix)
        h = string.gsub(h, '@name@', name)
        c = string.gsub(c1 .. c, '@PREFIX@', prefix)
        c = string.gsub(c, '@name@', name)
        stream_h = stream_h .. h
        stream_c = stream_c .. c
    end
    return stream_h, stream_c
end

function config.parse(schema, filename, runmain, globals, args)
-- Parse(run) a lua file
--
-- filename: the filename of the lua file
-- runmain: if true, the main function in the file is executed

    local namespace = setmetatable({}, {__index=globals})

    namespace['__file__'] = filename
    namespace['args'] = {}
    for i, k in pairs(args) do
        namespace['args'][i - 1] = k 
    end

    local param, err = loadfile(filename, 'bt', namespace)
    if param == nil then
        error(err)
    end
    param()

    -- prune main function from the parameter file
    local main = namespace['main']
    namespace['main'] = nil

    if main ~= nil then
        if runmain then
            main()
        end
    end

    local namespace2 = schema.bind(namespace)

    local ret, err = dump.tostring(namespace2)
    if ret == nil then
        error(err)
    end
    return ret

end

return config
