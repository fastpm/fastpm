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

typedef enum { FASTPM_MODEL_NONE, FASTPM_MODEL_LINEAR, FASTPM_MODEL_ZA, FASTPM_MODEL_2LPT, FASTPM_MODEL_PM } FastPMModelType;
typedef enum { FASTPM_KERNEL_3_4, FASTPM_KERNEL_3_2, FASTPM_KERNEL_5_4, FASTPM_KERNEL_EASTWOOD } FastPMKernelType;

typedef struct {
    /* input parameters */
    size_t nc;
    double boxsize;
    double omega_m;
    double alloc_factor;
    VPMInit * vpminit;
    int USE_COLA;
    int USE_DX1_ONLY;
    int USE_NONSTDDA;
    int USE_ZOLA;
    double nLPT;

    FastPMKernelType KERNEL_TYPE;
    FastPMModelType USE_MODEL;
    int K_LINEAR;

    /* Extensions */
    FastPMExtension * exts[12];

    /* internal variables */
    MPI_Comm comm;
    int ThisTask;
    int NTask;

    PMStore * p;
    VPM * vpm_list;

    FastPMModel * model;
    /* Pointer to the current PM object */
    PM * pm;
    PM * pm_2lpt;
} FastPM;

enum FastPMExtensionPoint {
    FASTPM_EXT_AFTER_FORCE,
    FASTPM_EXT_BEFORE_KICK,
    FASTPM_EXT_BEFORE_DRIFT,
    FASTPM_EXT_MAX,
};

typedef int 
    (* fastpm_ext_after_force) 
    (FastPM * fastpm, FastPMFloat * deltak, double a_x, void * userdata);
typedef int 
    (* fastpm_ext_before_kick) 
    (FastPM * fastpm, FastPMKick * kick, void * userdata);
typedef int
    (* fastpm_ext_before_drift) 
    (FastPM * fastpm, FastPMDrift * drift, void * userdata);

typedef struct FastPMExtension {
    void * function; /* The function signature must match the types above */
    void * userdata;
    struct FastPMExtension * next;
} FastPMExtension;

typedef struct FastPMDrift {
    FastPM * fastpm;
    PMStore * p;
    double dyyy;
    double da1;
    double da2;
    double af;
} FastPMDrift;

typedef struct FastPMKick {
    FastPM * fastpm;
    PMStore * p;
    float q1;
    float q2;
    float dda;
    double af;
} FastPMKick;

void fastpm_init(FastPM * fastpm, 
    int NprocY,  /* Use 0 for auto */
    int UseFFTW, /* Use 0 for PFFT 1 for FFTW */
    MPI_Comm comm);

void 
fastpm_add_extension(FastPM * fastpm, 
    enum FastPMExtensionPoint where,
    void * function, void * userdata);

void 
fastpm_destroy(FastPM * fastpm);

void 
fastpm_setup_ic(FastPM * fastpm, FastPMFloat * delta_k_ic);

void
fastpm_evolve(FastPM * fastpm, double * time_step, int nstep);

typedef int 
(*fastpm_interp_action) (FastPM * fastpm, PMStore * pout, double aout, void * userdata);

/* This function can be used in after_kick and after_drift plugins for 
 * interpolating and writing a snapshot */
void 
fastpm_interp(FastPM * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata);

void fastpm_drift_init(FastPMDrift * drift, FastPM * fastpm, PMStore * pi, double af);
void fastpm_kick_init(FastPMKick * kick, FastPM * fastpm, PMStore * pi, double af);
void
fastpm_kick_one(FastPMKick * kick, ptrdiff_t i, float vo[3]);
void
fastpm_drift_one(FastPMDrift * drift, ptrdiff_t i, double xo[3]);

FASTPM_END_DECLS
