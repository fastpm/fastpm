typedef struct FastPMModel FastPMModel;

struct FastPMModel {
    FastPMModelType type;
    FastPM * fastpm;
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
    void (*build)(FastPMModel * model, PMStore * p, double ainit, double afinal);
    void (*evolve)(FastPMModel * model, double af);
    void (*destroy)(FastPMModel * model);
};

double
fastpm_model_find_correction(FastPMModel * model,
    double a_x, double a_x1, double a_v, double a_v1);

void fastpm_model_init(FastPMModel * model, FastPM * fastpm, FastPMModelType type);
void fastpm_model_pt_init(FastPMModel * model);
void fastpm_model_pm_init(FastPMModel * model);
void fastpm_model_linear_init(FastPMModel * model);

void fastpm_model_create_subsample(FastPMModel * model, PMStore * psub, int attributes);
double fastpm_model_measure_large_scale_power(FastPMModel * model, PMStore * p);
void fastpm_model_destroy(FastPMModel * model);
void fastpm_model_build(FastPMModel * model, PMStore * p, double ainit, double afinal);
void fastpm_model_evolve(FastPMModel * model, double af);

void fastpm_calculate_forces(FastPM * fastpm, FastPMFloat * delta_k);

void 
fastpm_kick_store(FastPM * fastpm, 
              PMStore * pi, PMStore * po,
              double af);

void 
fastpm_drift_store(FastPM * fastpm,
               PMStore * pi, PMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPM * fastpm,
                PMStore * p, PMStore * po,
                double aout);
