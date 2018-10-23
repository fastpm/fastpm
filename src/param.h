typedef struct {
    int UseFFTW;
    int NprocY;
    int Nwriters;
    size_t MemoryPerRank;
    LuaConfig * config;
    char * string;
} Parameters;

void
free_parameters(Parameters * prr);

Parameters *
parse_args(int argc, char ** argv, char ** error);

#define CONF(prr, name) lua_config_get_ ## name (prr->config)
#define HAS(prr, name) lua_config_has_ ## name (prr->config)

#ifdef MPI_VERSION
Parameters *
parse_args_mpi(int argc, char ** argv, char ** error, MPI_Comm comm);
#endif


