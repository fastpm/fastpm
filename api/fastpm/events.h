enum FastPMEventType {
    FASTPM_EVENT_FORCE,
    FASTPM_EVENT_INTERPOLATION,
    FASTPM_EVENT_TRANSITION,
    FASTPM_EVENT_MAX,
};

enum FastPMEventStage {
    FASTPM_EVENT_STAGE_BEFORE,
    FASTPM_EVENT_STAGE_AFTER,
};

typedef struct {
    enum FastPMEventType type;
    enum FastPMEventStage stage;
} FastPMEvent;

typedef int
    (* FastPMEventHandlerFunction)(FastPMSolver * solver, FastPMEvent * event, void * userdata);

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
    FastPMFloat * delta_k;
    double a_f;
} FastPMForceEvent;

struct FastPMEventHandler {
    enum FastPMEventType type;
    enum FastPMEventStage stage;
    FastPMEventHandlerFunction function; /* The function signature must match the types above */
    void * userdata;
    struct FastPMEventHandler * next;
};

void 
fastpm_solver_add_event_handler(FastPMSolver * fastpm,
    enum FastPMEventType type,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function, void * userdata);

void
fastpm_solver_emit_event(FastPMSolver * fastpm,
    enum FastPMEventType where, enum FastPMEventStage stage, FastPMEvent * event);
