FASTPM_BEGIN_DECLS
typedef struct VPM VPM;

typedef struct VPMInit {
    double a_start;
    int pm_nc_factor;
} VPMInit;

typedef struct FastPMDrift FastPMDrift;
typedef struct FastPMKick FastPMKick;
typedef struct FastPMExtension FastPMExtension;
typedef struct FastPMModel FastPMModel;

typedef enum { FASTPM_FORCE_FASTPM = 0, FASTPM_FORCE_PM, FASTPM_FORCE_COLA} FastPMForceType;
typedef enum { FASTPM_MODEL_NONE, FASTPM_MODEL_LINEAR, FASTPM_MODEL_ZA, FASTPM_MODEL_2LPT, FASTPM_MODEL_PM } FastPMModelType;
typedef enum { FASTPM_KERNEL_3_4, FASTPM_KERNEL_3_2, FASTPM_KERNEL_5_4,
               FASTPM_KERNEL_GADGET,
               FASTPM_KERNEL_EASTWOOD,
               FASTPM_KERNEL_NAIVE,
            } FastPMKernelType;
typedef enum { FASTPM_DEALIASING_NONE,
               FASTPM_DEALIASING_GAUSSIAN, FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN,
               FASTPM_DEALIASING_TWO_THIRD, FASTPM_DEALIASING_GAUSSIAN36 } FastPMDealiasingType;

typedef struct {
    PM * pm;
    PMStore * p;

    MPI_Comm comm;
    int NTask;
    int ThisTask;

    /* input parameters */
    size_t nc;
    double boxsize;
    double omega_m;
    double alloc_factor;
    VPMInit * vpminit;
    int USE_DX1_ONLY;
    int USE_NONSTDDA;
    int USE_SHIFT;
    int SAVE_Q;

    double nLPT;
    FastPMPainterType PAINTER_TYPE;
    int painter_support;
    FastPMForceType FORCE_TYPE;
    FastPMKernelType KERNEL_TYPE;
    FastPMDealiasingType DEALIASING_TYPE;
    FastPMModelType USE_MODEL;
    int K_LINEAR;

    /* Extensions */
    FastPMExtension * exts[12];

    struct {
        /* For printing only. Do not use them to derive any physics quantities. */
        int istep;
        double a_x;
        double a_x1;
        double a_v;
        double a_v1;
        double dx1[3];
        double dx2[3];
        struct {
            double min;
            double max;
        } imbalance;
        int Nmesh;
    } info;

    VPM * vpm_list;

    FastPMModel * model;
    PM * basepm;
    FastPMPainter painter[1];
} FastPMSolver;

enum FastPMExtensionPoint {
    FASTPM_EXT_AFTER_FORCE,
    FASTPM_EXT_BEFORE_KICK,
    FASTPM_EXT_BEFORE_DRIFT,
    FASTPM_EXT_MAX,
};

typedef int 
    (* fastpm_ext_after_force) 
    (FastPMSolver * fastpm, FastPMFloat * deltak, double a_x, void * userdata);
typedef int 
    (* fastpm_ext_before_kick) 
    (FastPMSolver * fastpm, FastPMKick * kick, void * userdata);
typedef int
    (* fastpm_ext_before_drift) 
    (FastPMSolver * fastpm, FastPMDrift * drift, void * userdata);

struct FastPMExtension {
    void * function; /* The function signature must match the types above */
    void * userdata;
    struct FastPMExtension * next;
};

struct FastPMDrift {
    FastPMSolver * fastpm;
    PMStore * p;
    double dyyy;
    double da1;
    double da2;
    double af;
};

struct FastPMKick {
    FastPMSolver * fastpm;
    PMStore * p;
    float q1;
    float q2;
    float dda;
    double af;
};

void fastpm_init(FastPMSolver * fastpm, 
    int NprocY,  /* Use 0 for auto */
    int UseFFTW, /* Use 0 for PFFT 1 for FFTW */
    MPI_Comm comm);

void 
fastpm_add_extension(FastPMSolver * fastpm, 
    enum FastPMExtensionPoint where,
    void * function, void * userdata);

void 
fastpm_destroy(FastPMSolver * fastpm);

void 
fastpm_setup_ic(FastPMSolver * fastpm, FastPMFloat * delta_k_ic);

void
fastpm_evolve(FastPMSolver * fastpm, double * time_step, int nstep);

typedef int 
(*fastpm_interp_action) (FastPMSolver * fastpm, PMStore * pout, double aout, void * userdata);

/* This function can be used in after_kick and after_drift plugins for 
 * interpolating and writing a snapshot */
void 
fastpm_interp(FastPMSolver * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata);

void fastpm_drift_init(FastPMDrift * drift, FastPMSolver * fastpm, PMStore * pi, double af);
void fastpm_kick_init(FastPMKick * kick, FastPMSolver * fastpm, PMStore * pi, double af);
void
fastpm_kick_one(FastPMKick * kick, ptrdiff_t i, float vo[3]);
void
fastpm_drift_one(FastPMDrift * drift, ptrdiff_t i, double xo[3]);

double
fastpm_growth_factor(FastPMSolver * fastpm, double a);

FASTPM_END_DECLS
