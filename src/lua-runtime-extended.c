#include <unistd.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static int l_get_nprocs(lua_State *L)
{
    int flag = 0;
    MPI_Initialized(&flag);
    int np = 1;
    if(flag) {
        MPI_Comm_size(MPI_COMM_WORLD, &np);
    }
    lua_pushinteger(L, np);
    return 1;
}

static int l_get_nthreads(lua_State *L)
{
    int nt;
#ifdef _OPENMP
    nt = omp_get_max_threads();
#else
    nt = 1;
#endif
    lua_pushinteger(L, nt);
    return 1;
}

static int l_getcwd(lua_State *L) {
    char r[1024];
    if(NULL == getcwd(r, 1024)) {
        r[0] = '.';
        r[1] = 0;
    }
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

    /* get_nprocs */
    lua_pushstring(L, "get_nprocs");
    lua_pushcfunction(L, l_get_nprocs);
    lua_settable(L, -3);

    lua_pushstring(L, "get_nthreads");
    lua_pushcfunction(L, l_get_nthreads);
    lua_settable(L, -3);

    lua_pop(L, 1);
    return 0;
}
