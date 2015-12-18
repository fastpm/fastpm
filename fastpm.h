
typedef struct {
    MPI_Comm comm;
    int ThisTask;
    int NTask;

    int nc;
    double boxsize;
    double omega_m;
    PMStore * p;
    PM * pm_2lpt;
    double * time_step;
    int n_time_step;
    int USE_COLA;
    int USE_NONSTDDA;
    int USE_LINEAR_THEORY;
    double nLPT;
    double K_LINEAR;
    VPM * vpm_list;

    int istep;
    PM * pm;

} FastPM;

void 
fastpm_prepare_ic(FastPM * fastpm, FastPMFloat * delta_k);

void
fastpm_decompose(FastPM * fastpm);

void 
fastpm_kick(FastPM * fastpm, 
              PMStore * pi, PMStore * po,
              double af);

void 
fastpm_drift(FastPM * fastpm,
               PMStore * pi, PMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPM * fastpm,
                PMStore * p, PMStore * po,
                double aout);

double fastpm_growth_factor(FastPM * fastpm, double a);


typedef int 
(*fastpm_interp_action) (FastPM * fastpm, PMStore * pout, void * userdata);

void 
fastpm_interp(FastPM * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata);

void 
fastpm_set_time(FastPM * fastpm, 
    int istep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1);

int 
fastpm_read_runpb_ic(FastPM * fastpm, PMStore * p, char * filename);

int 
fastpm_write_runpb_snapshot(FastPM * fastpm, PMStore * p, char * filebase);
