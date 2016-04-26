# bash

produce_source () {
    xxd -i $1
}
produce_open () {
    local fn=$1
    local prefix=${fn//./_}
    local prefix=${prefix//-/_}
cat <<EOF
    {
        if(luaL_loadbuffer(L, (char*) $prefix, ${prefix}_len, "$fn")
        || lua_pcall(L, 0, 2, 0)
        ) {
            return 1;
        }
        const char * name = lua_tostring(L, -2);
        lua_setglobal(L, name);
        lua_pop(L, 1);
    }
EOF
}

# header
cat <<EOF
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

EOF
# definition of scripts
for script in $*; do
    produce_source $script
done


# entrance function
cat <<EOF

#include "lua-runtime-extended.c"

int lua_open_runtime(lua_State * L)
{
EOF
for script in $*; do
    produce_open $script
done
cat <<EOF
    lua_extend_runtime(L);
    return 0;
}
EOF
