#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>

#include <fastpm/logging.h>
#include <fastpm/string.h>


typedef struct FastPMMSGHandler FastPMMSGHandler;

struct FastPMMSGHandler {
    fastpm_msg_handler handler;
    void * userdata;
    MPI_Comm comm;
    FastPMMSGHandler * prev;
};

#if HAS_BACKTRACE
#include <execinfo.h>

static void fastpm_abort() {
    int size;
    void * ptrs[128];
    size = backtrace(&ptrs[0], 128);
    backtrace_symbols_fd(&ptrs[0], size, 0);
    abort();
}
#else
static void fastpm_abort() {
    abort();
}
#endif
static FastPMMSGHandler handler_data = {
    .handler = NULL,
    .userdata = NULL,
    .comm = (MPI_Comm) 0,
    .prev = NULL
};

void fastpm_void_msg_handler(
            const enum FastPMLogLevel level,
            const enum FastPMLogType type,
            const int errcode,
             const char * message, 
            MPI_Comm comm,
            void * userdata) {
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);
    if(type == COLLECTIVE) {
        MPI_Barrier(comm);
    }

    if(ThisTask != 0 && type == COLLECTIVE)
        return;

    if(level == ERROR) {
        if(type == COLLECTIVE) {
            if(ThisTask == 0) {
                fprintf(stdout, "%s", message);
                fflush(stdout);
            }
        } else {
            fprintf(stdout, "ThisTask = %d %s", ThisTask, message);
            fflush(stdout);
        }
        fastpm_abort();
    }
}

void fastpm_default_msg_handler(
            const enum FastPMLogLevel level,
            const enum FastPMLogType type,
            const int errcode,
             const char * message, 
            MPI_Comm comm,
            void * userdata) {
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);
    if(type == COLLECTIVE) {
        MPI_Barrier(comm);
    }

    if(type == COLLECTIVE) {
        if(ThisTask == 0) {
            fprintf(stdout, "%s", message);
            fflush(stdout);
        }
    } else {
        if(level == ERROR) {
            fprintf(stdout, "ThisTask = %d %s", ThisTask, message);
        } else {
            fprintf(stdout, "%s", message);
        }
        fflush(stdout);
    }

    if(ThisTask != 0 && type == COLLECTIVE)
        return;

    if(level == ERROR) fastpm_abort();
}

void fastpm_set_msg_handler(fastpm_msg_handler handler, MPI_Comm comm, void * userdata)
{
    handler_data.handler = handler;
    handler_data.userdata = userdata;
    handler_data.comm = comm;
}

void fastpm_push_msg_handler(fastpm_msg_handler handler, MPI_Comm comm, void * userdata)
{
    FastPMMSGHandler * prev = malloc(sizeof(FastPMMSGHandler));
    *prev = handler_data;
    handler_data.prev = prev;
    handler_data.handler = handler;
    handler_data.userdata = userdata;
    handler_data.comm = comm;
}

void fastpm_pop_msg_handler()
{
    FastPMMSGHandler * prev = handler_data.prev;
    handler_data = *prev;
    free(prev);
}

static double now()
{
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_sec + now.tv_usec * 1e-6;
}

static char * process(const char * buffer) {
    char header[128];

    static double t0 = -1.0;
    if(t0 < 0) t0 = now();
    sprintf(header, "[ %012.04f ]: ", now() - t0);

    int Nlines = 0;
    const char * p;
    for(p = buffer; *p; p ++) {
        if(*p == '\n') Nlines ++;
    }
    int need_newline = 0;
    /* at least 1 char and not ending with \n*/
    if(p != buffer && *(p - 1) != '\n') {
        need_newline = 1;
        Nlines ++;
    }
    char * ret = malloc(strlen(buffer) + strlen(header) * Nlines + need_newline + 1);
    char * q = ret;

    strcpy(q, header);
    q += strlen(header);

    for(p = buffer; *p; p ++) {
        *q = *p;
        q ++;
        if (*p == '\n' && *(p+1) != 0) {
            strcpy(q, header);
            q += strlen(header);
        }
    }
    if(need_newline) {
        *q = '\n';
        q ++;
    }
    *q = 0;
    return ret;
}


void fastpm_log2(const enum FastPMLogLevel level, 
            const enum FastPMLogType type,
            const int code, 
            const char * fmt, va_list argp) {

    if (handler_data.handler == NULL) {
        fastpm_set_msg_handler(fastpm_default_msg_handler, MPI_COMM_WORLD, NULL);
    }

    char * buffer = fastpm_strdup_vprintf(fmt, argp);
    char * processed = process(buffer);
    handler_data.handler(level, type, code, processed, handler_data.comm, handler_data.userdata);
    free(processed);
    free(buffer);
}

void fastpm_ilog(const enum FastPMLogLevel level, 
            const char * fmt, ...) {
    va_list argp;
    va_start(argp, fmt);
    fastpm_log2(level, INDIVIDUAL, 0, fmt, argp); 
    va_end(argp);
}

void fastpm_log(const enum FastPMLogLevel level, 
            const char * fmt, ...) {
    va_list argp;
    va_start(argp, fmt);
    fastpm_log2(level, COLLECTIVE, 0, fmt, argp); 
    va_end(argp);
}

void fastpm_info(const char * fmt, ...) {
    va_list argp;
    va_start(argp, fmt);
    fastpm_log2(INFO, COLLECTIVE, 0, fmt, argp); 
    va_end(argp);
}
void fastpm_raise(const int code, const char * fmt, ...) {
    va_list argp;
    va_start(argp, fmt);
    fastpm_log2(ERROR, INDIVIDUAL, code, fmt, argp); 
    va_end(argp);
}
