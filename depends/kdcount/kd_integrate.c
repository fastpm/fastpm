#include <alloca.h>
#include "kdtree.h"

typedef struct TraverseData {
    KDAttr * attr;
    uint64_t* count;
    double * weight;
    double * min;
    double * max;
    int attr_ndims;
    int node_ndims;
    uint64_t brute_force;
    uint64_t node_node;
} TraverseData;

static void
kd_integrate_check(TraverseData * trav, KDNode * node)
{

    ptrdiff_t i;
    int d;
    KDTree * t0 = node->tree;
    int node_ndims = trav->node_ndims;
    int attr_ndims = trav->attr_ndims;

    double * p0base = alloca(node->size * sizeof(double) * node_ndims);
    double * w0base = alloca(node->size * sizeof(double) * attr_ndims);
    /* collect all node positions to a continue block */
    double *p0, *w0;

    kd_collect(node, &t0->input, p0base);
    if(attr_ndims > 0) {
        kd_collect(node, &trav->attr->input, w0base);
    }

    for (p0 = p0base, w0 = w0base, i = 0; i < node->size; i++) {
        int inside = 1;
        for (d = 0; d < node_ndims; d++){
            if(p0[d] < trav->min[d] || p0[d] >= trav->max[d]) {
                inside = 0;
                break;
            }
        }
        if(inside) {
            trav->count[0] += 1;
            trav->brute_force ++;
            for(d = 0; d < attr_ndims; d++) {
                trav->weight[d] += w0[d];
            }
        }
        w0 += attr_ndims;
        p0 += node_ndims;
    }
}


static void 
kd_integrate_traverse(TraverseData * trav, KDNode * node) 
{
    int node_ndims = trav->node_ndims;
    int attr_ndims = trav->attr_ndims;
    int d;
    double *min0 = kd_node_min(node);
    double *max0 = kd_node_max(node);

    for(d = 0; d < node_ndims; d++) {
        if(min0[d] >= trav->max[d] || max0[d] < trav->min[d]) {
            /* fully outside, skip this node */
            return;
        }
    }

    int inside = 1;
    for(d = 0; d < node_ndims; d++) {
        if(min0[d] < trav->min[d] || max0[d] >= trav->max[d]) {
            inside = 0;
            break;
        }
    }
    if(inside) {
        /* node inside integration range */
        trav->count[0] += node->size;
        trav->node_node++;
        if(attr_ndims > 0) {
            double * w0 = kd_attr_get_node(trav->attr, node);
            for(d = 0; d < attr_ndims; d++) {
                trav->weight[d] += w0[d];
            }
        }
        return;
    }

    if(node->dim < 0) {
        /* can't open the node, need to enumerate */
        kd_integrate_check(trav, node);
    } else {
        kd_integrate_traverse(trav, node->link[0]);
        kd_integrate_traverse(trav, node->link[1]);
    } 
}

void 
kd_integrate(KDNode * node, KDAttr * attr,
        uint64_t * count, double * weight,
        double * min, double * max,
        uint64_t * brute_force,
        uint64_t * node_node)
{
    int attr_ndims;
    int d;

    if (attr)
        attr_ndims = attr->input.dims[1];
    else
        attr_ndims = 0;

    TraverseData trav = {
        .attr = attr,
        .count = count,
        .weight = weight,
        .attr_ndims = attr_ndims,
        .node_ndims = node->tree->input.dims[1],
        .min = min,
        .max = max,
        .brute_force = 0,
        .node_node = 0,
    };

    count[0] = 0;
    if(attr_ndims > 0) {
        for(d = 0; d < attr_ndims; d++) {
            weight[d] = 0;
        }
    }

    kd_integrate_traverse(&trav, node);

    *brute_force = trav.brute_force;
    *node_node = trav.node_node;

}
