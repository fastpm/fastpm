#include "kdtree.h"

struct TraverseData
{
    double maxr;
    double maxr2;
    double openr2;
    void * userdata;
    kd_enum_visit_edge visit_edge;
    kd_enum_check_nodes check_nodes;
    kd_enum_visit_node visit_node;
    void * check_nodes_data;
    int always_open;
    int skip_symmetric;
};

static int kd_enum_internal(struct TraverseData * trav, KDNode * nodes[2]);

static int _kd_enum_check_nodes(void * data, KDEnumNodePair * pair)
{
    struct TraverseData * trav = data;
    return kd_enum_check(pair->nodes, trav->maxr2, trav->skip_symmetric,
            trav->visit_edge, trav->userdata);
}

int
kd_enum(KDNode * nodes[2], double maxr,
        kd_enum_visit_edge visit_edge,
        kd_enum_check_nodes check_nodes,
        void * userdata)
{
    return kd_enum_full(nodes, maxr, visit_edge, check_nodes, NULL, 1.0, 0, userdata);
}

int
kd_enum_full(KDNode * nodes[2], double maxr,
        kd_enum_visit_edge visit_edge,
        kd_enum_check_nodes check_nodes,
        kd_enum_visit_node visit_node,
        double opening_factor,
        int skip_symmetric,
        void * userdata)
{
    struct TraverseData trav = {
        .maxr = maxr,
        .maxr2 = maxr * maxr,
        .check_nodes = check_nodes,
        .visit_edge = visit_edge,
        .userdata = userdata,
        .visit_node = visit_node,
        .openr2 = maxr * maxr * opening_factor * opening_factor,
        .skip_symmetric = skip_symmetric,
    };
    if(check_nodes == NULL) {
        trav.check_nodes = _kd_enum_check_nodes;
        trav.check_nodes_data = &trav;
    } else {
        trav.check_nodes_data = userdata;
    }
    return kd_enum_internal(&trav, nodes);

}
static int
_kd_enum_visit_node(struct TraverseData * trav, KDNode * node)
{
    if (trav->visit_node == NULL) return 1;
    return trav->visit_node(trav->userdata, node) != 0;
}
/*
 * enumerate two KDNode trees, up to radius max.
 *
 * for each pair i in nodes[0] and j in nodes[1],
 * if the distance is smaller than maxr,
 * call callback.
 * if callback returns nonzero, terminate and return the value
 * */
static int kd_enum_internal(struct TraverseData * trav, KDNode * nodes[2])
{
    int Nd = nodes[0]->tree->input.dims[1];
    double distmax = 0, distmin = 0;
    int d;
    double *min0 = kd_node_min(nodes[0]);
    double *min1 = kd_node_min(nodes[1]);
    double *max0 = kd_node_max(nodes[0]);
    double *max1 = kd_node_max(nodes[1]);
    for(d = 0; d < Nd; d++) {
        double min, max;
        double realmin, realmax;
        min = min0[d] - max1[d];
        max = max0[d] - min1[d];
        kd_realminmax(nodes[0]->tree, min, max, &realmin, &realmax, d);
        distmin += realmin * realmin;
        distmax += realmax * realmax;
    }
    /*
    printf("%g %g %g \n", distmin, distmax, maxr * maxr);
    print(nodes[0]);
    print(nodes[1]);
    */
    if (distmin > trav->maxr2 * 1.00001) {
        /* nodes are too far, skip them */
        return 0;
    }
    if (distmax >= trav->openr2) {
        /* nodes may intersect, open them */
        int open = nodes[0]->size < nodes[1]->size;
        if(nodes[open]->dim < 0 || !_kd_enum_visit_node(trav, nodes[open])) {
            open = (open == 0);
        }
        if(nodes[open]->dim < 0 || !_kd_enum_visit_node(trav, nodes[open])) {
            /* can't open the nodes, need to enumerate */
        } else {
            KDNode * save = nodes[open];
            nodes[open] = save->link[0];
            int rt;
            rt = kd_enum_internal(trav, nodes);
            if(rt != 0) {
                nodes[open] = save;
                return rt;
            }
            nodes[open] = save->link[1];
            rt = kd_enum_internal(trav, nodes);
            nodes[open] = save;
            return rt;
        }
    } else {
        /* fully inside, fall through,
         * and enumerate  */
    }

    KDEnumNodePair pair;
    pair.nodes[0] = nodes[0];
    pair.nodes[1] = nodes[1];
    pair.distmin2 = distmin;
    pair.distmax2 = distmax;

    return trav->check_nodes(trav->check_nodes_data, &pair);
}


int
kd_enum_check(KDNode * nodes[2], double maxr2, int skip_symmetric, kd_enum_visit_edge visit_edge, void * userdata)
{
    int rt = 0;
    ptrdiff_t i, j;
    int d;
    KDTree * t0 = nodes[0]->tree;
    KDTree * t1 = nodes[1]->tree;
    int Nd = t0->input.dims[1];

    double * p0base = malloc(nodes[0]->size * sizeof(double) * Nd);
    double * p1base = malloc(nodes[1]->size * sizeof(double) * Nd);
    /* collect all nodes[1] positions to a continue block */
    double * p1, * p0;

    KDEnumPair pair;

    /* no need to collect weight */
    kd_collect(nodes[0], &t0->input, p0base);
    kd_collect(nodes[1], &t1->input, p1base);

    for (p0 = p0base, i = nodes[0]->start; 
        i < nodes[0]->start + nodes[0]->size; i++, p0 += Nd) {
        pair.i = t0->ind[i];
        for (p1 = p1base, j = nodes[1]->start; 
             j < nodes[1]->start + nodes[1]->size; j++, p1 +=Nd) {
            pair.j = t1->ind[j];
            if (skip_symmetric && pair.i >= pair.j) continue;
            double r2 = 0.0;
            if (t0 != t1 || pair.i != pair.j) {
                for (d = 0; d < Nd; d++){
                    double dx = p1[d] - p0[d];
                    dx = kd_realdiff(nodes[0]->tree, dx, d);
                    r2 += dx * dx;
                }
            }
            if(r2 <= maxr2) {
                pair.r = sqrt(r2);
                if(0 != visit_edge(userdata, &pair)) {
                    rt = -1;
                    goto exit;
                }
            }
        }
    }
exit:
    free(p1base);
    free(p0base);
    return rt;
}

