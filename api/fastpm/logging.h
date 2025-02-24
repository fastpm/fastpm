#ifndef __FASTPM_LOGGING_H__
#define __FASTPM_LOGGING_H__

#ifndef FASTPM_BEGIN_DECLS
#ifdef __cplusplus
#define FASTPM_BEGIN_DECLS extern "C" {
#define FASTPM_END_DECLS }
#else
#define FASTPM_BEGIN_DECLS
#define FASTPM_END_DECLS
#endif
#endif

FASTPM_BEGIN_DECLS

enum FastPMLogLevel {
    ERROR = 100,
    INFO = 1,
};

enum FastPMLogType {
    COLLECTIVE = 0,
    INDIVIDUAL = 1,
};

typedef void 
(*fastpm_msg_handler)(
    const enum FastPMLogLevel level, 
    const enum FastPMLogType type, 
    const int errcode, 
    const char * message, 
    MPI_Comm comm,
    void * userdata);

void 
fastpm_ilog(const enum FastPMLogLevel level, 
            const char * fmt, ...);

void 
fastpm_log(const enum FastPMLogLevel level, 
           const char * fmt, ...);

void 
fastpm_set_msg_handler(fastpm_msg_handler handler, MPI_Comm comm, void * userdata);

void 
fastpm_push_msg_handler(fastpm_msg_handler handler, MPI_Comm comm, void * userdata);

void 
fastpm_pop_msg_handler();

void 
fastpm_default_msg_handler(
        const enum FastPMLogLevel level,
        const enum FastPMLogType type,
        const int errcode, 
        const char * message, 
        MPI_Comm comm,
        void * userdata);

void 
fastpm_void_msg_handler(
        const enum FastPMLogLevel level,
        const enum FastPMLogType type,
        const int errcode, 
        const char * message, 
        MPI_Comm comm,
        void * userdata);

void fastpm_info(const char * fmt, ...);
void fastpm_raise(const int code, const char * fmt, ...);

FASTPM_END_DECLS

#endif
