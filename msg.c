//
// Utility functions to write message, record computation time ...
//

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "msg.h"

static int myrank= -1;
static int Log_level;

// Initialize using msg_ functions.
void msg_init()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
}

void msg_set_loglevel(const enum LogLevel lv)
{
  Log_level= lv;
}

void msg_printf(const enum LogLevel msg_level, const char *fmt, ...)
{
  if(myrank == 0 && msg_level >= Log_level) {
    va_list argp;

    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp);
    fflush(stdout);
    va_end(argp);
  }
}

void msg_abort(const int errret, const char *fmt, ...)
{
  va_list argp;

  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  MPI_Abort(MPI_COMM_WORLD, errret);
}  

