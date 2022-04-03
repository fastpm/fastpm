#ifndef __FASTPM_HISTOGRAM_H__
#define __FASTPM_HISTOGRAM_H__

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

typedef struct {
    int64_t * counts;
    double * edges;
    int Nedges;
} FastPMHistogram;

void
fastpm_histogram_init(FastPMHistogram * hist, double amin, double amax, int nedges);

void 
fastpm_histogram_destroy(FastPMHistogram * hist);


FASTPM_END_DECLS

#endif
