#include "kdtree.h"

typedef struct TraverseData {
    ptrdiff_t * head;
    ptrdiff_t * next;
    ptrdiff_t * len;
    double ll;

    ptrdiff_t visited;
    ptrdiff_t connected;
    ptrdiff_t enumerated;
    ptrdiff_t maxdepth;
    ptrdiff_t nsplay;
    ptrdiff_t totaldepth;
} TraverseData;

static int 
_kd_fof_callback(void * data, KDEnumPair * pair) 
{
    TraverseData * trav = (TraverseData*) data;
    ptrdiff_t i = pair->i;
    ptrdiff_t j = pair->j;

    trav->visited++;

    if(trav->head[i] == trav->head[j]) return 0;

    if (trav->len[i] < trav->len[j] ) {
        /* merge in the shorter list */
        ptrdiff_t tmp;
        tmp = i;
        i = j;
        j = tmp;
    }

    /* update the length */
    trav->len[trav->head[i]] += trav->len[trav->head[j]];
    ptrdiff_t oldnext = trav->next[i];

    trav->next[i] = trav->head[j];

    ptrdiff_t k, kk;
    kk = trav->head[j]; /* shut-up the compiler, we know the loop will 
                          run at least once. */

    /* update the head marker of each element of the joined list */
    for(k = trav->head[j]; k >= 0 ; kk=k, k = trav->next[k]) {
        trav->head[k] = trav->head[i];
        trav->totaldepth++;
    }
    /* append items after i to the end of the merged list */
    trav->next[kk] = oldnext;
    //printf("reattaching %td to %td\n", oldnext, kk);

    return 0;
}

extern struct {
    /* performance counters */
    ptrdiff_t visited;
    ptrdiff_t enumerated;
    ptrdiff_t connected;
    ptrdiff_t maxdepth;
    ptrdiff_t nsplay;
    ptrdiff_t totaldepth;
} last_traverse;

int 
kd_fof_linkedlist(KDNode * node, double linking_length, ptrdiff_t * head)
{
    KDNode * nodes[2] = {node, node};
    TraverseData * trav = & (TraverseData) {};

    trav->head = head;
    trav->next = malloc(sizeof(trav->next[0]) * node->tree->input.dims[0]);
    trav->len = malloc(sizeof(trav->len[0]) * node->tree->input.dims[0]);
    trav->ll = linking_length;

    ptrdiff_t i;
    for(i = node->start; i < node->start + node->size; i ++) {
        ptrdiff_t j = node->tree->ind[i];
        trav->head[j] = j;
        trav->next[j] = -1;
        trav->len[j] = 1;
    }

    trav->visited = 0;
    trav->enumerated = 0;
    trav->maxdepth = 0;
    trav->totaldepth = 0;
    trav->nsplay = 0;
    trav->connected = 0;

    kd_enum_full(nodes, linking_length, _kd_fof_callback, NULL, NULL, 1.0, 1, trav);

    last_traverse.visited = trav->visited;
    last_traverse.enumerated = trav->enumerated;
    last_traverse.connected = trav->connected;
    last_traverse.maxdepth = trav->maxdepth;
    last_traverse.nsplay = trav->nsplay;
    last_traverse.totaldepth = trav->totaldepth;

    free(trav->next);
    free(trav->len);
    return 0;

exc_bad:
    free(trav->next);
    free(trav->len);
    return -1;
}

