FASTPM_BEGIN_DECLS
#define FASTPM_EVENT_FORCE "FORCE"
#define FASTPM_EVENT_LPT "LPT"
#define FASTPM_EVENT_TRANSITION "TRANSITION"
#define FASTPM_EVENT_INTERPOLATION "INTERPOLATION"

typedef struct VPM VPM;

typedef struct VPMInit {
    double a_start;
    int pm_nc_factor;
} VPMInit;

typedef struct FastPMDriftFactor FastPMDriftFactor;
typedef struct FastPMKickFactor FastPMKickFactor;

typedef struct {
    FastPMEvent base;
    FastPMDriftFactor * drift;
    FastPMKickFactor * kick;
    double a1;
    double a2;
} FastPMInterpolationEvent;

typedef struct {
    FastPMEvent base;
    FastPMTransition * transition;
} FastPMTransitionEvent;

typedef struct {
    FastPMEvent base;
    PM * pm;
    FastPMFloat * delta_k;
    FastPMStore * p;
} FastPMLPTEvent;

typedef struct {
    FastPMEvent base;
    FastPMKernelType kernel;
    FastPMPainter * painter;
    PM * pm;
    FastPMFloat * delta_k;
    double N; /* total number of particles painted. */
    double a_f;
    double a_n; /* time of next force calculation; or -1. if already the last force calculation. */
} FastPMForceEvent;

typedef struct {
    size_t nc;
    double boxsize;
    double omega_m;
    double hubble_param;
    double alloc_factor;

    VPMInit * vpminit;
    int USE_DX1_ONLY;
    int USE_SHIFT;

    FastPMColumnTags ExtraAttributes;

    double nLPT;
    /* FIXME: give them better looking names. */
    FastPMPainterType PAINTER_TYPE;
    int painter_support;
    FastPMForceType FORCE_TYPE;
    FastPMKernelType KERNEL_TYPE;
    FastPMDealiasingType DEALIASING_TYPE;

    int NprocY;  /* Use 0 for auto */
    int UseFFTW; /* Use 0 for PFFT 1 for FFTW */
} FastPMConfig;

#define FASTPM_SOLVER_NSPECIES 6

typedef struct {
    PM * pm;
    FastPMStore species[FASTPM_SOLVER_NSPECIES];
    char has_species[FASTPM_SOLVER_NSPECIES];     //interesting, use a character. kind of like an array of 0s and 1s saying if species exits....?
    //below keeps track of order species are added so that they can be freed in the correct order.
    char add_species_order[FASTPM_SOLVER_NSPECIES];  //element 0 is id of first species added, element 1 sis id of second species added...
    int N_added_species;    //the number of species currently added to the solver

    MPI_Comm comm;
    int NTask;
    int ThisTask;

    /* input parameters */
    FastPMConfig config[1];

    /* cosmology */
    FastPMCosmology cosmology[1];

    /* Extensions */
    FastPMEventHandler * event_handlers;

    VPM * vpm_list;

    PM * basepm;
} FastPMSolver;

enum FastPMAction {
    FASTPM_ACTION_FORCE,
    FASTPM_ACTION_KICK,
    FASTPM_ACTION_DRIFT,
};

struct FastPMDriftFactor {
    FastPMForceType forcemode;
    double ai;
    double ac;
    double af;

    /* */
    int nsamples;
    double Dv1; /* at ac */
    double Dv2; /* at ac */

    double dyyy[32];
    double da1[32];
    double da2[32];
};

struct FastPMKickFactor {
    FastPMForceType forcemode;

    double ai;
    double ac;
    double af;

    int nsamples;
    double q1;
    double q2;
    double dda[32];
    double Dv1[32];
    double Dv2[32];
};

void
fastpm_solver_init(FastPMSolver * fastpm, FastPMConfig * config, MPI_Comm comm);

void 
fastpm_solver_destroy(FastPMSolver * fastpm);

FastPMStore *
fastpm_solver_get_species(FastPMSolver * fastpm, enum FastPMSpecies species);

void
fastpm_solver_add_species(FastPMSolver * fastpm, enum FastPMSpecies species, size_t np_total);

void
fastpm_solver_destroy_species(FastPMSolver * fastpm);

void 
fastpm_solver_setup_lpt(FastPMSolver * fastpm, 
                enum FastPMSpecies species,
                FastPMFloat * delta_k_ic,
                double a0);

PM *
fastpm_find_pm(FastPMSolver * fastpm, double a);

void
fastpm_solver_evolve(FastPMSolver * fastpm, double * time_step, int nstep);

void fastpm_drift_init(FastPMDriftFactor * drift, FastPMSolver * fastpm, double ai, double ac, double af);
void fastpm_kick_init(FastPMKickFactor * kick, FastPMSolver * fastpm, double ai, double ac, double af);
void fastpm_kick_one(FastPMKickFactor * kick, FastPMStore * p,  ptrdiff_t i, float vo[3], double af);
void fastpm_drift_one(FastPMDriftFactor * drift, FastPMStore * p, ptrdiff_t i, double xo[3], double ae);

double
fastpm_solver_growth_factor(FastPMSolver * fastpm, double a);
double
fastpm_solver_growth_rate(FastPMSolver * fastpm, double a);

void 
fastpm_kick_store(FastPMKickFactor * kick,
    FastPMStore * pi, FastPMStore * po, double af);

void 
fastpm_drift_store(FastPMDriftFactor * drift,
               FastPMStore * pi, FastPMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPMSolver * fastpm,
                FastPMSolver * snapshot,
                FastPMDriftFactor * drift, FastPMKickFactor * kick,
                double aout);
void 
fastpm_set_species_snapshot(FastPMSolver * fastpm,
                FastPMStore * p,
                FastPMDriftFactor * drift, FastPMKickFactor * kick,
                FastPMStore * po,
                double aout);
void 
fastpm_unset_snapshot(FastPMSolver * fastpm,
                FastPMSolver * snapshot,
                FastPMDriftFactor * drift, FastPMKickFactor * kick,
                double aout);
void 
fastpm_unset_species_snapshot(FastPMSolver * fastpm,
                FastPMStore * p,
                FastPMDriftFactor * drift, FastPMKickFactor * kick,
                FastPMStore * po,
                double aout);

FASTPM_END_DECLS
