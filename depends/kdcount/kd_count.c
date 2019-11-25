#include <alloca.h>
#include "kdtree.h"

typedef struct TraverseData {
    KDAttr *attrs[2];
    int nedges;
    double * edges;
    uint64_t * count;
    double * weight;
    int node_ndims; 
    int attr_ndims; /* 0 if attrs[...] are NULL */
    kd_point_point_cullmetric ppcull;
    kd_node_node_cullmetric nncull;
    void * userdata;
    uint64_t brute_force;
    uint64_t node_node;
} TraverseData;

static inline int 
lower_bound(double key, double * r2, int N) 
{
    /* finding the first item that is equal or greater than key. */
    int left = 0, right = N;
    if(N == 0) return 0;
    if(key < r2[0]) return 0;
    if(key > r2[N-1]) return N;
    while(right > left) {
        int mid = left + ((right - left) >> 1);
        if(key > r2[mid]) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return left;
}
#if 0
This is not used.
static int bisect_right(double key, double * r2, int N) {
    if(key <= r2[0]) return 0;
    if(key > r2[N-1]) return N;
    int left = 0, right = N;
    while(right > left) {
        int mid = left + ((right - left) >> 1);
        if(key >= r2[mid]) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return left;
}
#endif
static void 
kd_count_check(TraverseData * trav, KDNode * nodes[2], 
        int start, int end) 
{

    ptrdiff_t i, j;
    int d;
    KDTree * t0 = nodes[0]->tree;
    KDTree * t1 = nodes[1]->tree;
    int node_ndims = trav->node_ndims;
    int attr_ndims = trav->attr_ndims;

    double * p0base = alloca(nodes[0]->size * sizeof(double) * node_ndims);
    double * p1base = alloca(nodes[1]->size * sizeof(double) * node_ndims);
    double * w0base = alloca(nodes[0]->size * sizeof(double) * attr_ndims);
    double * w1base = alloca(nodes[1]->size * sizeof(double) * attr_ndims);
    /* collect all nodes[1] positions to a continue block */
    double * p1, * p0, *w0, *w1;
    double * dx = alloca(sizeof(double) * node_ndims);

    kd_collect(nodes[0], &t0->input, p0base);
    kd_collect(nodes[1], &t1->input, p1base);
    if(attr_ndims > 0) {
        kd_collect(nodes[0], &trav->attrs[0]->input, w0base);
        kd_collect(nodes[1], &trav->attrs[1]->input, w1base);
    }
    for (p0 = p0base, w0 = w0base, i = 0; i < nodes[0]->size; i++) {
        for (p1 = p1base, w1 = w1base, j = 0; j < nodes[1]->size; j++) {
            double rr;
            for (d = 0; d < node_ndims; d++){
                dx[d] = kd_realdiff(nodes[0]->tree, p1[d] - p0[d], d);
            }
            int cull = trav->ppcull(trav->userdata, node_ndims, dx, &rr);
            int b = lower_bound(rr, &trav->edges[start], end - start) + start;
            if (!cull && b < trav->nedges) {
                trav->count[b] += 1;
                trav->brute_force ++;
                for(d = 0; d < attr_ndims; d++) {
                    trav->weight[b * attr_ndims + d] += w0[d] * w1[d];
                }
            }
            w1 += attr_ndims;
            p1 += node_ndims;
        }
        w0 += attr_ndims;
        p0 += node_ndims;
    }
    
}

static void 
kd_count_traverse(TraverseData * trav, KDNode * nodes[2], 
        int start, int end) 
{
    int node_ndims = trav->node_ndims;
    int attr_ndims = trav->attr_ndims;
    double distmax = 0, distmin = 0;
    int d;
    double *min0 = kd_node_min(nodes[0]);
    double *min1 = kd_node_min(nodes[1]);
    double *max0 = kd_node_max(nodes[0]);
    double *max1 = kd_node_max(nodes[1]);
    double *realmin = alloca(sizeof(double) * node_ndims);
    double *realmax = alloca(sizeof(double) * node_ndims);
    for(d = 0; d < node_ndims; d++) {
        realmin[d] = min0[d] - max1[d];
        realmax[d] = max0[d] - min1[d];
        kd_realminmax(nodes[0]->tree, realmin[d], realmax[d], &realmin[d], &realmax[d], d);
    }
    /* http://docs.cryengine.com/display/SDKDOC4/Culling+Explained */
    int cull = trav->nncull(trav->userdata, node_ndims, realmin, realmax, &distmin, &distmax);

    /* Two nodes are never going to make it to any binning */
    if(cull) return;

    start = lower_bound(distmin, &trav->edges[start], end - start) + start;
    end = lower_bound(distmax, &trav->edges[start], end - start) + start;
    if(start >= trav->nedges) {
        /* too far! skip */
        return;
    }
    if(start == end) {
        trav->node_node++;
        /* all bins are quickly counted no need to open*/
        trav->count[start] += nodes[0]->size * nodes[1]->size;
        if(attr_ndims > 0) {
            double * w0 = kd_attr_get_node(trav->attrs[0], nodes[0]);
            double * w1 = kd_attr_get_node(trav->attrs[1], nodes[1]);
            for(d = 0; d < attr_ndims; d++) {
                trav->weight[start * attr_ndims + d] += w0[d] * w1[d]; 
            }
        }
        return;
    }

    /* nodes may intersect, open them */
    int open = nodes[0]->size < nodes[1]->size;
    if(nodes[open]->dim < 0) {
        open = (open == 0);
    }
    if(nodes[open]->dim < 0) {
        /* can't open the nodes, need to enumerate */
        kd_count_check(trav, nodes, start, end);
    } else {
        KDNode * save = nodes[open];
        nodes[open] = save->link[0];
        kd_count_traverse(trav, nodes, start, end);
        nodes[open] = save->link[1];
        kd_count_traverse(trav, nodes, start, end);
        nodes[open] = save;
    } 
}

static int kd_node_node_euclidean(
    void * userdata,
    int ndims,
    double * realmin,
    double * realmax,
    double * distmin, double * distmax)
{
    int i;
    *distmin = 0;
    *distmax = 0;
    for(i = 0; i < ndims; i ++) {
        *distmin += realmin[i] * realmin[i];
        *distmax += realmax[i] * realmax[i];
    }
    return 0;
}

static int kd_point_point_euclidean(
    void * userdata,
    int ndims,
    double * dx,
    double * dist)
{
    int i;
    *dist = 0;
    for(i = 0; i < ndims; i ++) {
        *dist += dx[i] * dx[i];
    }
    return 0;
}

void
kd_count(KDNode * nodes[2], KDAttr * attrs[2], 
        double * edges, uint64_t * count, double * weight, 
        int nedges,
        kd_point_point_cullmetric ppcull,
        kd_node_node_cullmetric nncull,
        void * userdata,
        uint64_t * brute_force,
        uint64_t * node_node
) 
{
    double * edges2 = alloca(sizeof(double) * nedges);
    int attr_ndims;
    if (attrs[0]) 
        attr_ndims = attrs[0]->input.dims[1];
    else
        attr_ndims = 0;

    if (ppcull == NULL) {
        ppcull = kd_point_point_euclidean;
    }
    if (nncull == NULL) {
        nncull = kd_node_node_euclidean;
    }

    TraverseData trav = {
        .attrs = {attrs[0], attrs[1]},
        .nedges = nedges,
        .edges = edges2,
        .count = count,
        .weight = weight,
        .attr_ndims = attr_ndims,
        .node_ndims = nodes[0]->tree->input.dims[1],
        .ppcull = ppcull,
        .nncull = nncull,
        .userdata = userdata,
        .brute_force = 0,
        .node_node = 0,
    };

    int d;
    int i;
    for(i = 0; i < nedges; i ++) {
        if(edges[i] >= 0)
            edges2[i] = edges[i] * edges[i];
        else
            edges2[i] = i - nedges;
        count[i] = 0;
        if(attr_ndims > 0) {
            for(d = 0; d < attr_ndims; d++) {
                weight[i * attr_ndims + d] = 0;
            }
        }
    }
    
    kd_count_traverse(&trav, nodes, 0, nedges);
    *brute_force = trav.brute_force;
    *node_node = trav.node_node;
}
