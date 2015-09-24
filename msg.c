//
// Utility functions to write message, record computation time ...
//

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include "msg.h"

static int myrank= -1;
static int Log_level;
static double t0 = 0;
static double now()
{
    struct timeval tp;
    //struct timezone tzp;
    //gettimeofday(&tp,&tzp);
    gettimeofday(&tp, 0);

    return (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6;
}

// Initialize using msg_ functions.
void msg_init()
{
    t0 = now();
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
}

void msg_set_loglevel(const enum LogLevel lv)
{
    Log_level= lv;
}

static char * process(const char * fmt) {
    char buf[128];
    sprintf(buf, "[ %012.04f ]: ", now() - t0);
    char * ret = malloc(strlen(fmt) + strlen(buf) + 1);
    strcpy(ret, buf);
    strcat(ret, fmt);
    return ret;
}

void msg_printf(const enum LogLevel msg_level, const char *fmt, ...)
{
    if(myrank == 0 && msg_level >= Log_level) {
        va_list argp;
        va_start(argp, fmt);
        char * newfmt = process(fmt);
        vfprintf(stdout, newfmt, argp);
        free(newfmt);
        fflush(stdout);
        va_end(argp);
    }
}

void msg_aprintf(const enum LogLevel msg_level, const char *fmt, ...)
{
    if(msg_level >= Log_level) {
        va_list argp;

        va_start(argp, fmt);
        char * newfmt = process(fmt);
        vfprintf(stdout, newfmt, argp);
        free(newfmt);
        fflush(stdout);
        va_end(argp);
    }
}


void msg_abort(const int errret, const char *fmt, ...)
{
    va_list argp;

    va_start(argp, fmt);
    char * newfmt = process(fmt);
    vfprintf(stderr, fmt, argp);
    free(newfmt);
    va_end(argp);

    abort();
    MPI_Abort(MPI_COMM_WORLD, errret);
}  

