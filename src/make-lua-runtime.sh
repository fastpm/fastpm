# bash
if ! which xxd > /dev/null; then
    echo "No xxd found. install vim" >&2
    exit 1
fi
produce_source () {
    # if xxd is not found, die; most sane systems have xxd.
    xxd -i $1 || echo "You need the command xxd (installed with vim) to run this" >&2 && return 1
    return 0
}
produce_open () {
    local fn=$1
    local modulename=${fn//.lua}
    local prefix=${fn//./_}
    local prefix=${prefix//-/_}
cat <<EOF
    {
        luaL_getsubtable(L, LUA_REGISTRYINDEX, "_LOADED");
        if(luaL_loadbuffer(L, (char*) $prefix, ${prefix}_len, "$fn")) {
            return 1;
        }
        lua_pushstring(L, "$modulename");
        if(lua_pcall(L, 1, 1, 0)) {
            return 1;
        }
        lua_setfield(L, -2, "$modulename");
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
