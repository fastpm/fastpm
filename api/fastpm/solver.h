FASTPM_BEGIN_DECLS
typedef struct VPM VPM;

typedef struct VPMInit {
    double a_start;
    int pm_nc_factor;
} VPMInit;

typedef struct FastPMDriftFactor FastPMDriftFactor;
typedef struct FastPMKickFactor FastPMKickFactor;
typedef struct FastPMExtension FastPMExtension;

typedef struct {
    size_t nc;
    double boxsize;
    double omega_m;
    double alloc_factor;

    VPMInit * vpminit;
    int USE_DX1_ONLY;
    int USE_SHIFT;
    int SAVE_Q;
    int SAVE_POT;

    double nLPT;
    FastPMPainterType PAINTER_TYPE;
    int painter_support;
    FastPMForceType FORCE_TYPE;
    FastPMKernelType KERNEL_TYPE;
    FastPMDealiasingType DEALIASING_TYPE;
    int K_LINEAR;

    int NprocY;  /* Use 0 for auto */
    int UseFFTW; /* Use 0 for PFFT 1 for FFTW */
} FastPMConfig;

typedef struct {
    PM * pm;
    FastPMStore * p;

    MPI_Comm comm;
    int NTask;
    int ThisTask;

    /* input parameters */
    FastPMConfig config[1];

    /* gravity solver */
    FastPMGravity gravity[1];

    /* cosmology */
    FastPMCosmology cosmology[1];

    /* Extensions */
    FastPMExtension * exts[12];

    struct {
        /* For printing only. Do not use them to derive any physics quantities. */
        double dx1[3];
        double dx2[3];
        struct {
            double min;
            double max;
        } imbalance;
        int Nmesh;
    } info;

    VPM * vpm_list;

    PM * basepm;
} FastPMSolver;

enum FastPMExtensionPoint {
    FASTPM_EXT_AFTER_FORCE,
    FASTPM_EXT_INTERPOLATE,
    FASTPM_EXT_BEFORE_TRANSITION,
    FASTPM_EXT_MAX,
};

enum FastPMAction {
    FASTPM_ACTION_FORCE,
    FASTPM_ACTION_KICK,
    FASTPM_ACTION_DRIFT,
};

typedef int 
    (* fastpm_ext_after_force) 
    (FastPMSolver * fastpm, FastPMFloat * deltak, double a_x, void * userdata);

typedef int
    (* fastpm_ext_interpolate) 
    (FastPMSolver * fastpm, FastPMDriftFactor * drift, FastPMKickFactor * kick, double a1, double a2, void * userdata);

typedef int
    (* fastpm_ext_transition) 
    (FastPMSolver * fastpm, FastPMTransition * transition, void * userdata);

struct FastPMExtension {
    void * function; /* The function signature must match the types above */
    void * userdata;
    struct FastPMExtension * next;
};

struct FastPMDriftFactor {
    FastPMForceType forcemode;
    double ai;
    double ac;
    double af;

    /* */
    int nsamples;
    double dyyy[32];
    double da1[32];
    double da2[32];
    double Dv1; /* at ac */
    double Dv2; /* at ac */
};

struct FastPMKickFactor {
    FastPMForceType forcemode;

    double ai;
    double ac;
    double af;

    int nsamples;
    float q1;
    float q2;
    float dda[32];
    double Dv1[32];
    double Dv2[32];
};

void
fastpm_solver_init(FastPMSolver * fastpm, FastPMConfig * config, MPI_Comm comm);

void 
fastpm_solver_add_extension(FastPMSolver * fastpm, 
    enum FastPMExtensionPoint where,
    void * function, void * userdata);

void 
fastpm_solver_destroy(FastPMSolver * fastpm);

void 
fastpm_solver_setup_ic(FastPMSolver * fastpm, FastPMFloat * delta_k_ic);

void
fastpm_solver_evolve(FastPMSolver * fastpm, double * time_step, int nstep);

void fastpm_drift_init(FastPMDriftFactor * drift, FastPMSolver * fastpm, double ai, double ac, double af);
void fastpm_kick_init(FastPMKickFactor * kick, FastPMSolver * fastpm, double ai, double ac, double af);
void fastpm_kick_one(FastPMKickFactor * kick, FastPMStore * p,  ptrdiff_t i, float vo[3], double af);
void fastpm_drift_one(FastPMDriftFactor * drift, FastPMStore * p, ptrdiff_t i, double xo[3], double ae);

double
fastpm_solver_growth_factor(FastPMSolver * fastpm, double a);

void 
fastpm_kick_store(FastPMKickFactor * kick,
    FastPMStore * pi, FastPMStore * po, double af);

void 
fastpm_drift_store(FastPMDriftFactor * drift,
               FastPMStore * pi, FastPMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPMDriftFactor * drift, FastPMKickFactor * kick,
                FastPMStore * p, FastPMStore * po,
                double aout);

FASTPM_END_DECLS
