typedef struct {
    int UseFFTW;
    int NprocY;
    int Nwriters;
    int MaxThreads;
    char * RestartSnapshotPath;
    size_t MemoryPerRank;

    char ** argv;
    int argc;
} CLIParameters;

typedef struct{
    LuaConfig * config;
    char * string;
} LUAParameters;

void
free_lua_parameters(LUAParameters * prr);

void
free_cli_parameters(CLIParameters * prr);

CLIParameters *
parse_cli_args(int argc, char ** argv);

LUAParameters *
parse_config(char * filename, int argc, char ** argv, char ** error);

#define CONF(prr, name) lua_config_get_ ## name (prr->config)
#define HAS(prr, name) lua_config_has_ ## name (prr->config)

#ifdef MPI_VERSION
CLIParameters *
parse_cli_args_mpi(int argc, char ** argv, MPI_Comm comm);
LUAParameters *
parse_config_mpi(char * filename, int argc, char ** argv, char ** error, MPI_Comm comm);
#endif


