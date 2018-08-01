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
    fastpm_add_event_handler_free(handlers,
            where,
            stage,
            function,
            userdata, NULL);
}

void
fastpm_add_event_handler_free(FastPMEventHandler ** handlers,
    const char * where,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata, void (*free)(void * ptr))
{
    FastPMEventHandler * nh = malloc(sizeof(FastPMEventHandler));
    strncpy(nh->type, where, 32);
    nh->stage = stage;
    nh->function = function;
    nh->userdata = userdata;
    nh->free = free;
    nh->next = *handlers;
    *handlers = nh;
}

void
fastpm_remove_event_handler(FastPMEventHandler ** handlers,
    const char * where,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function,
    void * userdata)
{
    FastPMEventHandler * h0, * h1, * h2;

    for(h0 = NULL, h1 = *handlers;
        h1;
        h0 = h1, h1 = h2) {

        h2 = h1->next;

        if(h1->stage != stage) continue;
        if(0 != strcmp(h1->type, where)) continue;
        if(h1->function != function) continue;
        if(h1->userdata != userdata) continue;

        if(h1->free) {
            h1->free(h1->userdata);
        }

        free(h1);

        if(h0 == NULL) {
            /* remove head */
            *handlers = h2;
        } else {
            h0->next = h2;
        }
    }
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
    FastPMEventHandler * h, * h2;
    for(h = *handlers; h; h = h2) {
        h2 = h->next;
        if(h->free) {
            h->free(h->userdata);
        }
        free(h);
    }
    *handlers = NULL;
}
