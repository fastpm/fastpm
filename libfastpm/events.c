#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>

void 
fastpm_solver_add_event_handler(FastPMSolver * fastpm,
    enum FastPMEventType type,
    enum FastPMEventStage stage,
    FastPMEventHandlerFunction function, void * userdata)
{
    FastPMEventHandler * nh = malloc(sizeof(FastPMEventHandler));
    nh->type = type;
    nh->stage = stage;
    nh->function = function;
    nh->userdata = userdata;
    nh->next = fastpm->event_handlers;
    fastpm->event_handlers = nh;
}

void
fastpm_solver_emit_event(FastPMSolver * fastpm,
    enum FastPMEventType type, enum FastPMEventStage stage, FastPMEvent * event)
{
    /* fill in the mandatory items */
    event->type = type;
    event->stage = stage;
    FastPMEventHandler * handler = fastpm->event_handlers;
    for(; handler; handler = handler->next) {
        if(handler->type != type) continue;
        if(handler->stage != stage) continue;
        handler->function(fastpm, event, handler->userdata);
    }
}
