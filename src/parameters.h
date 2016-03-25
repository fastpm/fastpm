typedef struct {
    char * string;
    size_t nc;
    
    /* when to switch pm size, and to what */
    int * pm_nc_factor;
    double * change_pm;
    int n_pm_nc_factor;
    int n_change_pm;

    double np_alloc_factor;

    int random_seed;
    int inverted_ic;
    int remove_cosmic_variance;

    double boxsize;

    double omega_m, sigma8, h;

    double * aout; int n_aout;

    char * read_powerspectrum;
    char * write_powerspectrum;
    char * write_runpb_snapshot;
    char * write_snapshot;
    char * read_runpbic;
    char * read_grafic;
    char * read_noisek;
    char * read_noise;
    char * write_noisek;
    char * write_noise;

    int force_mode;
    int cola_stdda;
    int use_dx1_only;

    double * time_step; int n_time_step;
    int enforce_broadband;
    int enforce_broadband_mode;
    int enforce_broadband_kmax;
    int UseFFTW;
    int NprocY;
    int Nwriters;
} Parameters;

#define PM_MOND_NONE 0
#define PM_MOND_SIMPLE 1
#define PM_MOND_NBC 2

/* Important that modes with PM have the bit set */
#define FORCE_MODE_PM 8
#define FORCE_MODE_COLA 10

#define MODEL_LINEAR 1
#define MODEL_2LPT 2
