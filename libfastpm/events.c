#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>

/**
 * prepend an event handler to the linked list.
 *
 */
void
fastpm_add_event_handler(FastPMEventHandler ** handlers,
    const char * where,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata)
{
    FastPMEventHandler * nh = malloc(sizeof(FastPMEventHandler));
    strncpy(nh->type, where, 32);
    nh->stage = stage;
    nh->function = function;
    nh->userdata = userdata;
    nh->next = *handlers;
    *handlers = nh;
}

/**
 * emit an event on the list of handlers.
 *
 * the context is per event source
 * the userdata is per handler.
 */
void
fastpm_emit_event(FastPMEventHandler * handlers,
    const char * type, enum FastPMEventStage stage,
    FastPMEvent * event, void * context)
{
    /* fill in the mandatory items */
    strncpy(event->type, type, 31);
    event->stage = stage;
    FastPMEventHandler * handler = handlers;
    for(; handler; handler = handler->next) {
        if(0 != strcmp(handler->type, type)) continue;
        if(handler->stage != stage) continue;
        handler->function(context, event, handler->userdata);
    }
}

void
fastpm_destroy_event_handlers(FastPMEventHandler ** handlers)
{
    /* FIXME: free VPM and stuff. */
    FastPMEventHandler * h, * h2;
    for(h = *handlers; h; h = h2) {
        h2 = h->next;
        free(h);
    }
    *handlers = NULL;
}
