struct FastPMModel {
    FastPMModelType type;
    FastPMSolverPM * fastpm;
    FastPMPainter painter[1];
    PM * pm;
    double Pexpect;
    int factor;
    struct {
        PMStore * po;
        double a_x;
        double a_x1;
        double a_v;
        double a_v1;
    } ev;
    void * priv;
    void (*build)(FastPMModel * model, double ainit, double afinal);
    void (*evolve)(FastPMModel * model, double af);
    void (*destroy)(FastPMModel * model);
};

double
fastpm_model_find_correction(FastPMModel * model,
    double a_x, double a_x1, double a_v, double a_v1);

void fastpm_model_init(FastPMModel * model, FastPMSolverPM * fastpm, FastPMModelType type);
void fastpm_model_pt_init(FastPMModel * model);
void fastpm_model_pm_init(FastPMModel * model);
void fastpm_model_linear_init(FastPMModel * model);

void fastpm_model_create_subsample(FastPMModel * model, PMStore * psub);
double fastpm_model_measure_large_scale_power(FastPMModel * model, PMStore * p);
void fastpm_model_destroy(FastPMModel * model);
void fastpm_model_build(FastPMModel * model, double ainit, double afinal);
void fastpm_model_evolve(FastPMModel * model, double af);

void fastpm_calculate_forces(FastPMSolverPM * fastpm, FastPMFloat * delta_k);

void 
fastpm_kick_store(FastPMSolverPM * fastpm, 
              PMStore * pi, PMStore * po,
              double af);

void 
fastpm_drift_store(FastPMSolverPM * fastpm,
               PMStore * pi, PMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPMSolverPM * fastpm,
                PMStore * p, PMStore * po,
                double aout);
