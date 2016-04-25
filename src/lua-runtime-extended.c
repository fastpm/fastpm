#include <unistd.h>

static int l_getcwd(lua_State *L) {
    char r[1024];
    getcwd(r, 1024);
    lua_pushstring(L, r);
    return 1;
}

int lua_extend_runtime(lua_State * L)
{
    lua_getglobal(L, "os");
    /* getcwd */
    lua_pushstring(L, "getcwd");
    lua_pushcfunction(L, l_getcwd);
    lua_settable(L, -3);

    lua_pop(L, 1);
    return 0;
}
