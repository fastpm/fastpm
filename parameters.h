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

    double a_init;
    double a_final;
    double boxsize;

    double omega_m, sigma8, h;

    double * zout; int n_zout;

    char * power_spectrum_filename;
    char * measure_power_spectrum_filename; 
    char * snapshot_filename;
    char * readic_filename;

    int force_mode;
    int cola_stdda;
    double smoothing;
    int diff_order;
    int time_step;
} Parameters;

#define FORCE_MODE_ZA 0
#define FORCE_MODE_2LPT 1
#define FORCE_MODE_PM 2
#define FORCE_MODE_COLA 3

#define TIME_STEP_LOGA 0
#define TIME_STEP_A 1
#define TIME_STEP_GROWTH 2

int read_parameters(const int argc, char * argv[], 
        Parameters * param);

#endif
