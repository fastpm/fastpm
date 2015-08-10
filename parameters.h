#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef struct {
    size_t nc;
    int pm_nc_factor1;
    int pm_nc_factor2;
    float change_pm;

    double np_alloc_factor;

    int random_seed;
    int nrealization;

    int loglevel;

    double boxsize;

    double omega_m, sigma8, h;

    double * zout; int n_zout;

    int pm_mond_mode;
    double * pm_mond_parameters; int n_pm_mond_parameters;

    char * power_spectrum_filename;
    char * measure_power_spectrum_filename; 
    char * snapshot_filename;
    char * readic_filename;

    int force_mode;
    int cola_stdda;
    double smoothing;
    int diff_order;
    int poisson_order;

    double * time_step; int n_time_step;
    int enforce_broadband;
} Parameters;

#define PM_MOND_NONE 0
#define PM_MOND_SIMPLE 1
#define PM_MOND_NBC 2

#define FORCE_MODE_ZA 1
#define FORCE_MODE_2LPT 2
/* Important that modes with PM have the bit set */
#define FORCE_MODE_PM 8
#define FORCE_MODE_COLA1 9
#define FORCE_MODE_COLA 10

#define TIME_STEP_LOGA 0
#define TIME_STEP_A 1
#define TIME_STEP_GROWTH 2

int read_parameters(const int argc, char * argv[], 
        Parameters * param);

#endif
