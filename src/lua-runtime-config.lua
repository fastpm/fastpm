config = {}

local function Schema()
    local self = {}

    local names = {}
    local trace_actions = {}

    function self.declare(options)
        local name = options.name
        local type = options.type
        local required = options.required or false
        local default = options.default
        local choices = options.choices
        local action = options.action

        local subtype = options.type:match('array:(%a+)')
        if subtype ~= nil then
            type = 'array'
        end
        self[name] = {type=type, subtype=subtype, required=required, default=default, choices=choice, action=action}
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
        for i, v in pairs(value) do
            if type(v) ~= entry.subtype then
                error(string.format("entry %d of key `%s` is of type `%s`, but `%s` is expected",
                    i, name, type(v), entry.subtype))
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
                return entry.default
            end
        else

            check_type(name, entry, value)

            if entry.action ~= nil then
                entry.action(value)
            end
            return value
        end
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
                luaL_eval(lc->L, "@name@");
                const char * val = lua_tostring(lc->L, -1);
                lua_pop(lc->L, 1);
                if(val) return _strdup(val);
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
        ]])
    else
        return process([[
        @CTYPE@ * @PREFIX@_get_@name@_full(LuaConfig * lc, int * size)
        {
            luaL_eval(lc->L, "@name@");
            const int n = luaL_len(lc->L, -1);
            double * array = (@CTYPE@*) malloc(sizeof(@CTYPE@) * n);
            int i;
            for(i = 1; i <= n; ++i) {
                lua_pushinteger(lc->L, i);
                lua_gettable(lc->L, -2);
                @CTYPE@ x = @CONVERTOR@(lc->L, -1);
                lua_pop(lc->L,1);
                array[i-1] = x;
            }
            lua_pop(lc->L, 1);
            *size = n;
            return array;
        }
        @CTYPE@ * @PREFIX@_get_@name@(LuaConfig * lc)
        {
            int size;
            return @PREFIX@_get_@name@_full(lc, &size);
        }
        int @PREFIX@_get_n_@name@(LuaConfig * lc)
        {
            int size;
            void * a = @PREFIX@_get_@name@_full(lc, &size);
            free(a);
            return size;
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

local function compile(schema, opt)
    local preample_h = [[
        typedef struct LuaConfig LuaConfig;

        void lua_config_free(LuaConfig * lc);
        LuaConfig * lua_config_new(const char * luastring);
        const char * lua_config_error(LuaConfig * lc);
    ]]

    local preample_c = [[
    #include <lua.h>
    #include <lauxlib.h>
    #include <lualib.h>
    #include <string.h>
    #include <stdlib.h>

    typedef struct LuaConfig LuaConfig;

    struct LuaConfig {
        lua_State * L;
        char * error;
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
        if(luaL_eval(L, luastring)) {
            lc->error = _strdup(lua_tostring(L, -1));
        }
        return lc;
    }

    void lua_config_free(LuaConfig * lc)
    {
        if(lc->error) free(lc->error);
        lua_close(lc->L);
        free(lc);
    }
    const char * lua_config_error(LuaConfig * lc)
    {
        return lc->error;
    }

    ]]

    local stream_h = rtrim(preample_h)
    local stream_c = rtrim(preample_c)

    for i, fn in pairs(opt.global_headers) do
        local include = "#include <" .. fn .. ">\n"
        stream_c = stream_c .. include
    end
    for i, fn in pairs(opt.local_headers) do
        local include = "#include \"" .. fn .. "\"\n"
        stream_c = stream_c .. include
    end

    for name, entry in schema.pairs() do
        local h = rtrim(visit(name, entry, 'h'))
        local c = rtrim(visit(name, entry, 'c'))

        h = string.gsub(h, '@PREFIX@', opt.prefix)
        h = string.gsub(h, '@name@', name)
        c = string.gsub(c, '@PREFIX@', opt.prefix)
        c = string.gsub(c, '@name@', name)
        stream_h = stream_h .. h
        stream_c = stream_c .. c
    end
    return stream_h, stream_c
end

config.compile = compile

return 'config', config
