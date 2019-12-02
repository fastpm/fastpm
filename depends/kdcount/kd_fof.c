#include <alloca.h>
#include "kdtree.h"

/* Friend of Friend:
 *
 * Connected component via edge enumeration.
 *
 * The connected components are stored as trees.
 *
 * (visit) two vertices i, j connected by an edge.
 *         if splay(i) differ from splay(j), the components shall be
 *         visited, call merge(i, j)
 *
 * (merge) Join two trees if they are connected by adding the root
 *         of i as a child of root of j.
 *
 * (splay) Move a leaf node to the direct child of the tree root
 *         and returns the root.
 *
 * One can show this algorithm ensures splay(i) is an label of
 * max connected components in the graph.
 *
 * Suitable for application where finding edges of a vertice is more expensive
 * than enumerating over edges.
 *
 * In KD-Tree implementation we apply optimization for overdense regions
 *
 * Before FOF we connect all nodes that are less than linking length in size (internally connected)
 * 
 * (connect) connect all points in a node into a tree, the first particle is the root.
 * 
 * In FOF, we use the dual tree algorithm for edge enumeration, but if two nodes are minimaly sperated
 * by linking length and both are internally connected
 *
 *   (nodemerge) merge(node1->first, node2->first)
 * 
 * For typically low resolution cosmological simulations this optimization improves speed by a few percent.
 * The improvement increases for highly clustered data.
 *
 * The storage is O(N) for the output labels + O(M) for the connection status of kdtree nodes.
 *
 * */

typedef struct TraverseData {
    ptrdiff_t * head;
    ptrdiff_t * ind; /* tree->ind */
    char * node_connected;
    double ll;
    double ll2;

    /* performance counters */
    ptrdiff_t visited;
    ptrdiff_t connected;
    ptrdiff_t enumerated;
    ptrdiff_t maxdepth;
    ptrdiff_t nsplay;
    ptrdiff_t totaldepth;
    int safe;
    int allpairs;
    int heuristics;
    int prefernodes;
    int buggy;
} TraverseData;

typedef struct{
    TraverseData * trav;
    int self_connected;
    ptrdiff_t visited;
} VisitEdgeData;

static ptrdiff_t splay(TraverseData * d, ptrdiff_t i)
{
    ptrdiff_t depth = 0;
    ptrdiff_t r = i;
    /* First find the root */
    while(d->head[r] != r) {
        depth ++;
        r = d->head[r];
    }

    if(d->safe) {
        /* safe guard */
        while(d->head[i] != i) {
            ptrdiff_t t = d->head[i];
            d->head[i] = r;
            i = t;
        }
    } else {
        /* link the nodes directly to the root to keep the tree flat */
        d->head[i] = r;
    }

    /* update performance counters */
    if(depth > d->maxdepth) {
        d->maxdepth = depth;
    }
    d->totaldepth += depth;
    d->nsplay ++;

    return r;
}

static int
_kd_fof_visit_edge(void * data, KDEnumPair * pair);

static int
_kd_fof_visit_edge_connected(void * data, KDEnumPair * pair);

static int
_kd_enum_check_connected(KDNode * nodes[2], double maxr2, int skip_symmetric, kd_enum_visit_edge visit_edge, void * userdata);

static int
_kd_fof_check_nodes(void * data, KDEnumNodePair * pair)
{
    TraverseData * trav = (TraverseData*) data;

    if(trav->node_connected[pair->nodes[0]->index]
    && trav->node_connected[pair->nodes[1]->index]
    && !trav->allpairs
    )
    {
        trav->connected += (pair->nodes[0]->size * pair->nodes[1]->size); /* this count is duplicated. shall divide by two */
        if(trav->buggy) {
            /*this was the older buggy implemention that uses arbitary pair. */
            KDEnumPair epair;
            epair.r = sqrt(pair->distmin2);
            epair.i = trav->ind[pair->nodes[0]->start];
            epair.j = trav->ind[pair->nodes[1]->start];

            _kd_fof_visit_edge(data, &epair);
            trav->enumerated += 1;
            return 0;
        }
        /* two fully connected nodes are linked, simply link the first particle.  */
        if(trav->heuristics) {
            /* the enum function will count enumeration */
            _kd_enum_check_connected(pair->nodes, trav->ll2, 1, _kd_fof_visit_edge, data);
            return 0;
        } else {
            kd_enum_check(pair->nodes, trav->ll2, 1, _kd_fof_visit_edge_connected, data);
            trav->enumerated += (pair->nodes[0]->size * pair->nodes[1]->size); /* this count is duplicated. shall divide by two */
            return 0;
        }
    } else {
        kd_enum_check(pair->nodes, trav->ll2, 1, _kd_fof_visit_edge, data);
        trav->enumerated += (pair->nodes[0]->size * pair->nodes[1]->size); /* this count is duplicated. shall divide by two */
    }

    return 0;
}

static int
_kd_fof_visit_edge(void * data, KDEnumPair * pair) 
{
    TraverseData * trav = (TraverseData *) data;

    ptrdiff_t i = pair->i;
    ptrdiff_t j = pair->j;

    trav->visited ++;

    ptrdiff_t root_i = splay(trav, i);
    ptrdiff_t root_j = splay(trav, j);

    /* merge root_j as direct subtree of the root */
    /* this is also correct if root_j == root_i */
    trav->head[root_i] = root_j;

    /* terminate immediately if two nodes are self-connected and
     * we have linked a pair*/
    return 0;
}

static int
_kd_fof_visit_edge_connected(void * data, KDEnumPair * pair) 
{
    _kd_fof_visit_edge(data, pair);
    return -1;
}

static double kd_node_maxdist2(KDNode * node)
{
    int d;
    double dist = 0;
    for(d = 0; d < node->tree->input.dims[1]; d ++) {
        double dx = kd_node_max(node)[d] - kd_node_min(node)[d];
        dist += dx * dx;
    }
    return dist;
}

static void
connect(TraverseData * trav, KDNode * node, int parent_connected)
{
    int c = 0;
    if(parent_connected) {
        c = 1;
    } else {
        if(kd_node_maxdist2(node) <= trav->ll2) {
            ptrdiff_t i;
            ptrdiff_t r = trav->ind[node->start];
            for(i = node->start + 1; i < node->size + node->start; i ++) {
                trav->head[trav->ind[i]] = r;
            }
            c = 1;
        }
    }
    trav->node_connected[node->index] = c;

    if(node->dim != -1) {
        connect(trav, node->link[0], c);
        connect(trav, node->link[1], c);
    }
}
static int
_kd_fof_visit_node(void * data, KDNode * node)
{
    TraverseData * trav = (TraverseData*) data;
    if (trav->prefernodes) return 1;
    if (trav->allpairs) return 1;
    if (trav->buggy) return 1;
    /* if the nodes are very nearby and all connected, then there is no need to split the node*/
    return trav->node_connected[node->index] == 0;
}

struct {
    /* performance counters */
    ptrdiff_t visited;
    ptrdiff_t enumerated;
    ptrdiff_t connected;
    ptrdiff_t maxdepth;
    ptrdiff_t nsplay;
    ptrdiff_t totaldepth;
} last_traverse = {0};

static int 
kd_fof_internal(KDNode * node, double linking_length, ptrdiff_t * head, int safe, int allpairs, int heuristics, int prefernodes, int buggy)
{
    KDNode * nodes[2] = {node, node};
    TraverseData * trav = & (TraverseData) {};

    trav->head = head;
    trav->ll = linking_length;
    trav->ll2 = linking_length * linking_length;
    trav->node_connected = calloc(node->tree->size, 1);
    trav->ind = node->tree->ind;
    trav->safe = safe;
    trav->allpairs = allpairs;
    trav->prefernodes = prefernodes;
    trav->heuristics = heuristics;
    trav->buggy = buggy;
    ptrdiff_t i;
    for(i = node->start; i < node->start + node->size; i ++) {
        ptrdiff_t j = trav->ind[i];
        trav->head[j] = j;
    }

    trav->visited = 0;
    trav->enumerated = 0;
    trav->maxdepth = 0;
    trav->totaldepth = 0;
    trav->nsplay = 0;
    trav->connected = 0;

    connect(trav, node, 0);

    kd_enum_full(nodes, linking_length, NULL, _kd_fof_check_nodes, _kd_fof_visit_node, 1.0, 1, trav);

    for(i = node->start; i < node->start + node->size; i ++) {
        ptrdiff_t j = trav->ind[i];
        trav->head[j] = splay(trav, j);
    }

    free(trav->node_connected);

    last_traverse.visited = trav->visited;
    last_traverse.enumerated = trav->enumerated;
    last_traverse.connected = trav->connected;
    last_traverse.maxdepth = trav->maxdepth;
    last_traverse.nsplay = trav->nsplay;
    last_traverse.totaldepth = trav->totaldepth;
    return 0;
}

int 
kd_fof(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 1, 0, 0, 0, 0);
}

int
kd_fof_allpairs(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 1, 1, 0, 0, 0);
}

int
kd_fof_prefernodes(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 1, 0, 0, 1, 0);
}

int
kd_fof_unsafe(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 0, 0, 0, 0, 0);
}

int
kd_fof_heuristics(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 1, 0, 1, 0, 0);
}

int
kd_fof_buggy(KDNode * node, double linking_length, ptrdiff_t * head)
{
    return kd_fof_internal(node, linking_length, head, 1, 0, 0, 0, 1);
}


void
kd_fof_get_last_traverse_info(ptrdiff_t *visited, ptrdiff_t *enumerated, ptrdiff_t *connected,
                              ptrdiff_t *maxdepth, ptrdiff_t *nsplay, ptrdiff_t *totaldepth)
{
    *visited = last_traverse.visited;
    *enumerated = last_traverse.enumerated;
    *connected = last_traverse.connected;
    *maxdepth = last_traverse.maxdepth;
    *nsplay = last_traverse.nsplay;
    *totaldepth = last_traverse.totaldepth;
}

/* check and stop after visiting the first edge */

/* This function is not used by default because although in cases it reduces the number of enumeration by half,
 * the speed up is only a few percent. Apparently visit is so much more expensive than distance calculation*/

static ptrdiff_t
_kd_fof_guess_closest(KDNode * node, double * p0, ptrdiff_t size, int Nd)
{
    double * min = kd_node_min(node);
    double * max = kd_node_max(node);
    double * center = alloca(sizeof(double) * Nd);
    
    ptrdiff_t i;
    int d;
    for (d = 0; d < Nd; d ++) {
        center[d] = 0.5 * (min[d] + max[d]);
    }
    ptrdiff_t i_min = -1;
    double r2_min = 0;
    for(i = 0; i < size; i ++, p0 += Nd) {
        double r2 = 0;
        for(d = 0; d < Nd; d++) {
            double dx = p0[d] - center[d];
            r2 += dx * dx;
        }
        if (r2 < r2_min || i_min == -1) {
            i_min = i;
            r2_min = r2;
        }
    }
    return i_min;
}
static int
_kd_enum_check_connected(KDNode * nodes[2], double maxr2, int skip_symmetric, kd_enum_visit_edge visit_edge, void * userdata)
{
    TraverseData * trav = (TraverseData*) userdata;

    ptrdiff_t i, j;
    int d;
    KDTree * t0 = nodes[0]->tree;
    KDTree * t1 = nodes[1]->tree;
    int Nd = t0->input.dims[1];

    double * p0base = malloc(nodes[0]->size * sizeof(double) * Nd);
    double * p1base = malloc(nodes[1]->size * sizeof(double) * Nd);
    /* collect all nodes[1] positions to a continue block */
    double * p1, * p0;

    if(nodes[0]->size == 0 || nodes[1]->size == 0) return 0;

    KDEnumPair pair;

    if (nodes[0] == nodes[1]) {
        /* node with itself, do nothing */
        goto exit;
    }

    kd_collect(nodes[0], &t0->input, p0base);
    kd_collect(nodes[1], &t1->input, p1base);

    i = _kd_fof_guess_closest(nodes[1], p0base, nodes[0]->size, Nd);
    j = _kd_fof_guess_closest(nodes[0], p1base, nodes[1]->size, Nd);

    trav->enumerated += nodes[0]->size;
    trav->enumerated += nodes[1]->size;

    pair.i = t0->ind[nodes[0]->start + i];
    pair.j = t1->ind[nodes[1]->start + j];

    double r2 = 0.0;
    for (d = 0; d < Nd; d++){
        double dx = p1base[Nd * j + d] - p0base[Nd * i + d];
        dx = kd_realdiff(nodes[0]->tree, dx, d);
        r2 += dx * dx;
    }
    if(r2 <= maxr2) {
        pair.r = sqrt(r2);
        visit_edge(userdata, &pair);
        goto exit;
    }

    for (p0 = p0base, i = nodes[0]->start; 
        i < nodes[0]->start + nodes[0]->size; i++, p0 += Nd) {
        pair.i = t0->ind[i];
        for (p1 = p1base, j = nodes[1]->start; 
             j < nodes[1]->start + nodes[1]->size; j++, p1 +=Nd) {
            pair.j = t1->ind[j];
            trav->enumerated ++;
            if (pair.i >= pair.j) continue;
            double r2 = 0.0;
            for (d = 0; d < Nd; d++){
                double dx = p1[d] - p0[d];
                dx = kd_realdiff(nodes[0]->tree, dx, d);
                r2 += dx * dx;
            }
            if(r2 <= maxr2) {
                pair.r = sqrt(r2);
                visit_edge(userdata, &pair);
                goto exit;
            }
        }
    }

exit:

    free(p1base);
    free(p0base);
    return 0;
}
