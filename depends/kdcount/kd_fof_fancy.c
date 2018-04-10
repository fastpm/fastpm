#include "kdtree.h"
/*
 *
 * First identify all nodes node_max ** 2 < ll ** 2,
 *     Set particles' delegate to the node
 *
 * Dual tree walk:
 *     If node pairs mindist < ll ** 2
 *         Merge two nodes
 *     If node pairs maxdist > ll ** 2
 *         Skip
 *     Open non-leaf nodes
 *
 * If both nodes are leaf:
 *     For particle pairs < ll ** 2
 *         Merge two particles
 *
 * Merge two particles:
 *     Find delegate
 * */
typedef struct TraverseData {
    int node_ndims; 
    double linkinglength;
} TraverseData;

void 
kd_fof(KDNode * node, double linkinglength, ptrdiff_t * head) 
{
    TraverseData trav = {
        .node_ndims = nodes[0]->tree->input.dims[1],
        .linkinglength = linkinglength,
    };

    kd_fof_traverse(&trav, node);
}

static void 
kd_fof_traverse(TraverseData * trav, KDNode * nodes[2]) 
{
    int node_ndims = trav->node_ndims;
    double distmax = 0, distmin = 0;
    int d;
    double *min0 = kd_node_min(nodes[0]);
    double *min1 = kd_node_min(nodes[1]);
    double *max0 = kd_node_max(nodes[0]);
    double *max1 = kd_node_max(nodes[1]);
    for(d = 0; d < node_ndims; d++) {
        double min, max;
        double realmin, realmax;
        min = min0[d] - max1[d];
        max = max0[d] - min1[d];
        kd_realminmax(nodes[0]->tree, min, max, &realmin, &realmax, d);
        distmin += realmin * realmin;
        distmax += realmax * realmax;
    }

    if(distmax <= trav->linkinglength2) {
        update_nodes_head(trav, nodes);
        return;
    }
    int open = nodes[0]->size < nodes[1]->size;
    if(nodes[open]->dim < 0) {
        open = (open == 0);
    }
    if(nodes[open]->dim < 0) {
        /* can't open the nodes, need to enumerate */
        kd_fof_check(trav, nodes);
    } else {
        KDNode * save = nodes[open];
        nodes[open] = save->link[0];
        kd_fof_traverse(trav, nodes);
        nodes[open] = save->link[1];
        kd_fof_traverse(trav, nodes);
        nodes[open] = save;
    } 

}

static void 
kd_fof_check(TraverseData * trav, KDNode * nodes[2]) 
{
    ptrdiff_t i, j;
    int d;
    KDTree * t0 = nodes[0]->tree;
    KDTree * t1 = nodes[1]->tree;
    int node_ndims = trav->node_ndims;

    double * p0base = alloca(nodes[0]->size * sizeof(double) * node_ndims);
    double * p1base = alloca(nodes[1]->size * sizeof(double) * node_ndims);
    /* collect all nodes[1] positions to a continue block */
    double * p1, * p0;
    double half[node_ndims];
    double full[node_ndims];

    if(t0->boxsize) {
        for(d = 0; d < node_ndims; d++) {
            half[d] = t0->boxsize[d] * 0.5;
            full[d] = t0->boxsize[d];
        }
    }

    kd_collect(nodes[0], &t0->input, p0base);
    kd_collect(nodes[1], &t1->input, p1base);
    if(attr_ndims > 0) {
        kd_collect(nodes[0], &trav->attrs[0]->input, w0base);
        kd_collect(nodes[1], &trav->attrs[1]->input, w1base);
    }
    for (p0 = p0base, w0 = w0base, i = 0; i < nodes[0]->size; i++) {
        for (p1 = p1base, w1 = w1base, j = 0; j < nodes[1]->size; j++) {
            double rr = 0.0;
            for (d = 0; d < node_ndims; d++){
                double dx = p1[d] - p0[d];
                dx = kd_realdiff(nodes[0]->tree, dx, d);
                rr += dx * dx;
            }
            if(rr < trav->linkinglength2) {
                /* */
                head0 = get_head(i + nodes[0].start);
                head1 = get_head(j + nodes[1].start);
                merge(head0, head1);
            }
            w1 += attr_ndims;
            p1 += node_ndims;
        }
        w0 += attr_ndims;
        p0 += node_ndims;
    }
    

}


