#include "kdtree.h"

typedef struct TraverseData {
    KDAttr * xmass;
    KDAttr * mass;
    double * pos;
    double r_cut2;
    double eta2;
    int node_ndims;
    double * force;
    kd_force_func func;
    void * userdata;
    int64_t node_probed;
    int64_t node_computed;
    int64_t pair_computed;
} TraverseData;

/* compute x - y*/
static double
distance(KDTree * tree, double * x, double * y, double * dx)
{
    double r2 = 0;
    int Nd = tree->input.dims[1];
    int d;
    double half;
    for(d = 0; d < Nd; d++) {
        dx[d] = y[d] - x[d];
        if (tree->boxsize) {
            half = 0.5 * tree->boxsize[d];
            if (dx[d] > half) dx[d] = dx[d] - tree->boxsize[d];
            else if (dx[d] < -half) dx[d] = dx[d] + tree->boxsize[d];
        }
        r2 += dx[d] * dx[d];
    }
    return r2;
}

