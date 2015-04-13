#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef struct {
  size_t nc;
  int pm_nc_factor1;
  int pm_nc_factor2;
  float change_pm;

  double np_alloc_factor;

  int ntimestep;
  int random_seed;
  int nrealization;

  int loglevel;

  double a_final;
  double boxsize;

  double omega_m, sigma8, h;

  double* zout; int n_zout;

  /*
  char power_spectrum_file[64];
  char fof_filename[64];
  char snapshot_filename[64];
  char subsample_filename[64];
  char cgrid_filename[64];
  char init_filename[64];
  */
  char* power_spectrum_filename; int strlen_power_spectrum_filename;
  char* measure_power_spectrum_filename; int strlen_measure_power_spectrum_filename;
  char* fof_filename;        int strlen_fof_filename;
  char* snapshot_filename;   int strlen_snapshot_filename;
  char* subsample_filename;  int strlen_subsample_filename;
  char* cgrid_filename;      int strlen_cgrid_filename;
  char* init_filename;       int strlen_init_filename;
  char* readic_filename;   int strlen_readic_filename;
  // 'strlen' includes the null character. It is the number of chars.


  //int subsample_factor;
  double fof_linking_factor;
  double subsample_factor;
  int cgrid_nc;

  int write_longid; 
  int qpm;
  int nopm;
  int stdda;
  double smoothing;
  int diff_order;
  int time_step;
} Parameters;

#define TIME_STEP_LOGA 0
#define TIME_STEP_A 1
#define TIME_STEP_GROWTH 2
//
//int read_parameters(const char filename[], Parameters* const param);
int read_parameters(const int argc, char * argv[], 
		    Parameters* const param);


#endif
