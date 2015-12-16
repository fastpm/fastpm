#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef struct {
    size_t nc;
    
    /* when to switch pm size, and to what */
    int * pm_nc_factor;
    double * change_pm;
    int n_pm_nc_factor;
    int n_change_pm;

    double np_alloc_factor;

    int random_seed;
    int nrealization;

    int loglevel;

    double boxsize;

    double omega_m, sigma8, h;

    double * zout; int n_zout;

    char * power_spectrum_filename;
    char * measure_power_spectrum_filename; 
    char * snapshot_filename;
    char * readic_filename;
    char * readnoise_filename;

    int force_mode;
    int cola_stdda;

    double * time_step; int n_time_step;
    int enforce_broadband; 
    double enforce_broadband_kmax;
    int UseFFTW;
    int NprocY;
} Parameters;

#define PM_MOND_NONE 0
#define PM_MOND_SIMPLE 1
#define PM_MOND_NBC 2

/* Important that modes with PM have the bit set */
#define FORCE_MODE_PM 8
#define FORCE_MODE_COLA 10

#endif
