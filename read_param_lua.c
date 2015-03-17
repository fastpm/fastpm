#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>
#include "parameters.h"

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "msg.h"

static int myrank_;

int read_parameter_file(const char filename[], Parameters* const param);
static void bcast_string(char** string, int* len);
static void bcast_array_double(double** parray, int* len);

//int read_parameters(const char filename[], Parameters* const param)
int read_parameters(const int argc, char * argv[],
		    Parameters* const param)
{
  if(argc < 2)
    msg_abort(1, "Error: Parameter file not specified. cola_code param.lua\n");

  char const * const filename= argv[argc-1];

  //int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
  if(myrank_ == 0) {
    int ret= read_parameter_file(filename, param);
    if(ret != 0)
      msg_abort(1001, "Error: Unable to read parameter file: %s\n", filename);
  }

  // Share parameters with other nodes
  MPI_Bcast(param, sizeof(Parameters), MPI_BYTE, 0, MPI_COMM_WORLD);

  bcast_string(&param->power_spectrum_filename, 
	       &param->strlen_power_spectrum_filename);
  bcast_string(&param->measure_power_spectrum_filename, 
	       &param->strlen_measure_power_spectrum_filename);
  bcast_string(&param->fof_filename,        &param->strlen_fof_filename);
  bcast_string(&param->snapshot_filename,   &param->strlen_snapshot_filename);
  bcast_string(&param->subsample_filename,  &param->strlen_subsample_filename);
  bcast_string(&param->cgrid_filename,      &param->strlen_cgrid_filename);
  bcast_string(&param->init_filename,       &param->strlen_init_filename);

  bcast_array_double(&param->zout,          &param->n_zout);

  return 0;
}

static int read_int(lua_State* L, const char name[])
{
  lua_getglobal(L, name);
  if(!lua_isnumber(L, -1)) {
    msg_abort(1010, "Error: Parameter %s not found in the parameter file\n",
	      name);
  }

  int n = lua_tointeger(L, -1);
  lua_pop(L, 1);

  return n;
}

static double read_double(lua_State* L, const char name[])
{
  lua_getglobal(L, name);
  if(!lua_isnumber(L, -1)) {
    msg_abort(1020, "Error: parameter %s not found in the parameter file\n",
	    name);
  }

  double val = lua_tonumber(L, -1);
  lua_pop(L, 1);

  return val;
}

/*
static void read_string_fixed(lua_State* L, const char name[], 
			 char* val, const int n)
{
  lua_getglobal(L, name);
  if(!lua_isstring(L, -1)) {
    msg_printf(debug, "Parameter %s not found, possibly on optional one.\n", name);
    val[0]= '\0';
    return;
  }
  
  char const * const str= lua_tostring(L, -1);

  strncpy(val, str, n);

  lua_pop(L, 1);
}
*/

static char* read_string2(lua_State* L, const char name[], int* len, 
			  bool required)
{
  lua_getglobal(L, name);

  if(!lua_isstring(L, -1)) { // string name not found
    if(required)
      msg_abort(1030, "Error: Parameter %s not found in the parameter file\n",
		name);
    else {
      *len= 0;
      return 0;
    }
  }
  
  char const * const str= lua_tostring(L, -1);
  const int n= strlen(str) + 1; // +1 for the terminating null '\0'
  char* const val= malloc(sizeof(char)*n); assert(val);

  strncpy(val, str, n);

  lua_pop(L, 1);

  *len= n;
  return val;
}

  

static int read_bool(lua_State* L, const char name[])
{
  lua_getglobal(L, name);
  if(!lua_isboolean(L, -1)) {
    msg_abort(1030, "Error: Parameter %s not found in the parameter file\n",
	      name);
  }

  int n = lua_toboolean(L, -1);
  lua_pop(L, 1);

  return n;
}

static double* read_array_double(lua_State* L, const char name[], int *len)
{
  lua_getglobal(L, name);
  if(!lua_istable(L, -1)) {
    msg_abort(1031, "Error: Parameter %s not found or not an array in the parameter file\n");
  }

  //n= lua_objlen(L, -1); // old version
  const int n= luaL_len(L, -1);     // version 5.2-work3 (5.2.2_1 now) and later

  double* const array= (double*) malloc(sizeof(double)*n); assert(array);

  for(int i=1; i<=n; ++i) {
    lua_pushinteger(L, i);
    lua_gettable(L, -2);
    double x= lua_tonumber(L, -1);
    lua_pop(L,1);

    array[i-1]= x;
  }
  lua_pop(L, 1); // pop table "aout"

  *len= n;
  return array;
}

int read_parameter_file(const char filename[], Parameters* const param)
{
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  if(luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)) {
    fprintf(stderr, "%s", lua_tostring(L, -1));
    return -1;
  }

  param->nc= read_int(L, "nc");
  param->boxsize= read_double(L, "boxsize");

  param->a_final= read_double(L, "a_final");
  param->ntimestep= read_int(L, "ntimestep");
  param->zout= read_array_double(L, "output_redshifts", &param->n_zout);

  param->random_seed= read_int(L, "random_seed");
  param->nrealization= read_int(L, "nrealization");

  param->omega_m= read_double(L, "omega_m");
  param->h= read_double(L, "h");
  param->sigma8= read_double(L, "sigma8");

  param->pm_nc_factor= read_int(L, "pm_nc_factor");
  param->np_alloc_factor= read_double(L, "np_alloc_factor");
  param->loglevel= read_int(L, "loglevel");

  // File Names and optional parameters realated
  param->power_spectrum_filename=
    read_string2(L, "powerspectrum", &param->strlen_power_spectrum_filename, 
		 true);

  param->measure_power_spectrum_filename=
    read_string2(L, "measure_power", 
        &param->strlen_measure_power_spectrum_filename, false);

  param->fof_filename= 
    read_string2(L, "fof", &param->strlen_fof_filename, false);

  if(param->fof_filename) {
    param->fof_linking_factor= read_double(L, "linking_factor");
    // TODO: loglevel is not set yet, so, don't write too much now
    // echo parameter data later...
    msg_printf(verbose, "fof linking_factor read: %f\n", 
	       param->fof_linking_factor);
  }
  else
    param->fof_linking_factor= 0.0;


  param->snapshot_filename=
    read_string2(L, "snapshot", &param->strlen_snapshot_filename, false);
  
  param->subsample_filename=
    read_string2(L, "subsample", &param->strlen_subsample_filename, false);

  if(param->subsample_filename) {
    param->subsample_factor= read_double(L, "subsample_factor");
    msg_printf(verbose, "subsample_factor read: %d\n", param->subsample_factor);
  }
  else
    param->subsample_factor= 0;

  
  param->cgrid_filename=
    read_string2(L, "coarse_grid", &param->strlen_cgrid_filename, false);

  if(param->cgrid_filename) {
    param->cgrid_nc= read_int(L, "coarse_grid_nc");
    msg_printf(verbose, "coarse_grid_nc read: %d\n", param->cgrid_nc);
  }
  else
    param->cgrid_nc= 0;

  param->init_filename=
    read_string2(L, "initial", &param->strlen_init_filename, false);

  if(param->snapshot_filename || param->init_filename) {
    param->write_longid= read_bool(L, "write_longid");
    msg_printf(verbose, "write_longid= %d\n", param->write_longid);
  }

    /* true to run QPM instead of COLA. */
  param->qpm = read_bool(L, "qpm");
    /* 0 to turn off smoothing; in cells */
  param->smoothing = read_double(L, "smoothing");
  param->diff_order = read_int(L, "diff_order");
  param->loga_step= read_bool(L, "loga_step");


  lua_close(L);

  return 0;
}

void bcast_string(char** pstring, int* len)
{
  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert(ret1 == MPI_SUCCESS);

  const int n= *len;

  if(n == 0) {
    *pstring= 0;
    return;
  }

  if(myrank_ != 0) {
    *pstring= malloc(sizeof(char)*n);
  }
    assert(*pstring);

  const int ret2= MPI_Bcast(*pstring, n, MPI_CHAR, 0, MPI_COMM_WORLD);
    assert(ret2 == MPI_SUCCESS);
}

void bcast_array_double(double** parray, int* len)
{
  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert(ret1 == MPI_SUCCESS);

  const int n= *len;

  if(n == 0) {
    *parray= 0;
    return;
  }

  if(myrank_ != 0) {
    *parray= malloc(sizeof(double)*n);
  }
    assert(*parray);

  const int ret2= MPI_Bcast(*parray, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    assert(ret2 == MPI_SUCCESS);
}
