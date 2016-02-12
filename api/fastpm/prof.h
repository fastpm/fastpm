#ifndef __FASTPM_PROF_H__
#define __FASTPM_PROF_H__
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

typedef struct FastPMClock FastPMClock;

void fastpm_clock_in(FastPMClock * clock);
void fastpm_clock_out(FastPMClock * clock);

FastPMClock * 
fastpm_clock_find(const char * file, const char * func, const char * name);
void fastpm_clock_out_barrier(FastPMClock * clock, MPI_Comm comm);

#define CLOCK(name) FastPMClock * CLK ## name = fastpm_clock_find(__FILE__, __func__, # name);\
                    fastpm_clock_in(CLK ## name);

#define ENTER(name) fastpm_clock_in(CLK ## name)
#define LEAVE(name) fastpm_clock_out(CLK ## name)
#define LEAVEB(name, comm) fastpm_clock_out_barrier(CLK ## name, comm)

void fastpm_clock_stat(MPI_Comm comm);
void fastpm_report_memory(MPI_Comm comm);

FASTPM_END_DECLS
#endif
