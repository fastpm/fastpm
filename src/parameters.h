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

    double omega_m, h;

    double * aout; int n_aout;

    char * read_powerspectrum;

    double sigma8;
    double scalar_amp;
    double pivot_scalar;
    double scalar_spectral_index;
    double f_nl;
    char * f_nl_type;

    char * write_powerspectrum;
    char * write_nonlineark;
    char * write_runpb_snapshot;
    char * write_snapshot;
    char * read_runpbic;
    char * read_grafic;

    char * read_lineark;
    char * write_lineark;
    char * read_whitenoisek;
    char * write_whitenoisek;

    char * force_mode;
    char * kernel_type;
    char * dealiasing_type;

    int cola_stdda;
    int use_zola;
    int use_dx1_only;

    double * time_step; int n_time_step;

    char * enforce_broadband_mode;
    int enforce_broadband_kmax;
    int UseFFTW;
    int NprocY;
    int Nwriters;
} Parameters;
