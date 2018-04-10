typedef struct FastPMEventHandler FastPMEventHandler;

enum FastPMEventStage {
    FASTPM_EVENT_STAGE_BEFORE,
    FASTPM_EVENT_STAGE_AFTER,
};

typedef struct {
    char type[32];
    enum FastPMEventStage stage;
} FastPMEvent;

typedef int
    (* FastPMEventHandlerFunction)(void * context, FastPMEvent * event, void * userdata);

struct FastPMEventHandler {
    char type[32];
    enum FastPMEventStage stage;
    FastPMEventHandlerFunction function; /* The function signature must match the types above */
    void * userdata;
    struct FastPMEventHandler * next;
    void (*free) (void*);
};

void
fastpm_add_event_handler(FastPMEventHandler ** handlers,
    const char * type,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata);

/** add a event handler that frees the userdata upon removal **/
void
fastpm_add_event_handler_free(FastPMEventHandler ** handlers,
    const char * where,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata, void (*free)(void * ptr));

void
fastpm_remove_event_handler(FastPMEventHandler ** handlers,
    const char * type,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata);

void
fastpm_emit_event(FastPMEventHandler * handlers,
    const char * type, enum FastPMEventStage stage,
    FastPMEvent * event, void * context);

void
fastpm_destroy_event_handlers(FastPMEventHandler ** handlers);
