
typedef struct FastPMModel {
    FastPMModelType type;
    FastPM * fastpm;
    PM * pm;
    PMStore * psub;
    double a_x;
    double Pexpect;
    int factor;
    struct {
        PMStore * po;
        double a_x;
        double a_x1;
        double a_v;
        double a_v1;
    } ev;
} FastPMModel;

double
fastpm_model_find_correction(FastPMModel * model, 
    double a_x, double a_x1, double a_v, double a_v1);

void fastpm_model_init(FastPMModel * model, FastPM * fastpm, FastPMModelType type);
void fastpm_model_destroy(FastPMModel * model);
void fastpm_model_build(FastPMModel * model, PMStore * p, double ainit);
void fastpm_model_evolve(FastPMModel * model, double af);
double fastpm_growth_factor(FastPM * fastpm, double a);


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
