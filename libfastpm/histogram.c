#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/histogram.h>

void
fastpm_histogram_init(FastPMHistogram * hist, double amin, double amax, int nedges)
{
    int i;
    double * edges = malloc(sizeof(double) * nedges);

    for(i = 0; i < nedges - 1; i ++) {
        edges[i] = (amax - amin) * i / (nedges - 1) + amin;
    }

    edges[nedges - 1] = amax;
    hist->edges = edges;
    hist->Nedges = nedges;
    hist->counts = calloc(hist->Nedges + 1, sizeof(int64_t));
}

void 
fastpm_histogram_destroy(FastPMHistogram * hist)
{
    free(hist->counts);
    free(hist->edges);
}

