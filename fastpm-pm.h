typedef struct {
    /* input parameters */
    int nc;
    double boxsize;
    double omega_m;
    double * time_step;
    int n_time_step;

    int USE_COLA;
    int USE_NONSTDDA;
    int USE_LINEAR_THEORY;
    double nLPT;
    double K_LINEAR;

    /* internal variables */
    MPI_Comm comm;
    int ThisTask;
    int NTask;

    PMStore * p;
    PM * pm_2lpt;
    VPM * vpm_list;

    int istep;
    PM * pm;

} FastPM;

void fastpm_init(FastPM * fastpm, 
    double alloc_factor, 
    int NprocY, 
    int UseFFTW, 
    int * pm_nc_factor, 
    int n_pm_nc_factor, 
    double * change_pm,
    MPI_Comm comm);

void 
fastpm_destroy(FastPM * fastpm);
void 
fastpm_prepare_ic(FastPM * fastpm, FastPMFloat * delta_k);

typedef int (* fastpm_after_kick_action) (FastPM * fastpm, void * userdata);
typedef int (* fastpm_after_drift_action) (FastPM * fastpm, void * userdata);
typedef int (* fastpm_after_force_action) (FastPM * fastpm, FastPMFloat * deltak, double a_x, void * userdata);

void
fastpm_evolve(FastPM * fastpm, 
    fastpm_after_force_action after_force,
    fastpm_after_drift_action after_drift,
    fastpm_after_kick_action after_kick,
    void * userdata);

typedef int 
(*fastpm_interp_action) (FastPM * fastpm, PMStore * pout, double aout, void * userdata);

void 
fastpm_interp(FastPM * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata);

