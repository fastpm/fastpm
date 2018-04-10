#include "kdtree.h"

static void * 
kd_malloc(KDTree * tree, size_t size) 
{
    if(tree->malloc != NULL) {
        return tree->malloc(tree->userdata, size);
    } else {
        return malloc(size);
    }
}

static KDNode * 
kd_alloc(KDTree * tree) 
{
    KDNode * ptr = kd_malloc(tree, sizeof(KDNode) + 
            sizeof(double) * 2 * tree->input.dims[1]
            );
    ptr->link[0] = NULL;
    ptr->link[1] = NULL;
    ptr->tree = tree;
    return ptr;
}
static void
kd_build_update_min_max(KDNode * node, double min[], double max[])
{
    int Nd = node->tree->input.dims[1];
    int d;
    ptrdiff_t i;

    /* no leaf, directly compute from input */
    if(node->size > 0) {
        for(d = 0; d < Nd; d++) {
            max[d] = kd_input(node->tree, node->start + 0, d);
            min[d] = max[d];
        }
    } else {
        for(d = 0; d < Nd; d++) {
            max[d] = 0;
            min[d] = 0;
        }
    }
    for(i = 0; i < node->size; i++) {
        for (d = 0; d < Nd; d++) {
            double x = kd_input(node->tree, node->start + i, d);
            if (max[d] < x) max[d] = x;
            if (min[d] > x) min[d] = x;
        }
    }

}
static ptrdiff_t
kd_build_pivot(KDNode * node, int dim, double split)
{
    ptrdiff_t p, q;
    p = node->start;
    q = node->start + node->size - 1;
    while(p <= q) {
        if(kd_input(node->tree, p, dim) < split) {
            p ++;
        } else if(kd_input(node->tree, q, dim) >= split) {
            q --;
        } else {
            kd_swap(node->tree, p, q); 
            p ++;
            q --;
        }
    }
    /* invariance: data[<p] < split and data[>q] >= split
     * after loop p > q.
     * thus data[0...,  p-1] < split
     * data[q + 1...., end] >= split
     * p - 1 < q + 1
     * p < q + 2
     * p > q
     * thus p = q + 1 after the loop.
     *
     * 0 -> p -1 goes to left
     * and p -> end goes to right, so that
     *  left < split
     *  and right >= split
     *  */
    return p;
}
static int
kd_build_find_longest_dim(int Nd, double * min, double * max)
{
    double longest = max[0] - min[0];
    int dim = 0;
    int d;
    for(d = 1; d < Nd; d++) {
        double tmp = max[d] - min[d];
        if(tmp > longest) {
            dim = d;
            longest = tmp;
        }
    }
    return dim;
}
static int
kd_build_is_poor_split(KDNode * node, ptrdiff_t p)
{
    return (p == node->start || p == node->start + node->size);
}

static ptrdiff_t
kd_build_split(KDNode * node, double minhint[], double maxhint[], ptrdiff_t next) 
{
    KDTree * tree = node->tree;
    int d;
    int Nd = tree->input.dims[1];
    double * max = kd_node_max(node);
    double * min = kd_node_min(node);

    if(node->size <= tree->thresh) {
        /* do not split */
        kd_build_update_min_max(node, min, max);
        return next;
    }

    /* initialize the min max to the hints. They will be updated to the actual value in the end. */
    for(d = 0; d < Nd; d++) {
        max[d] = maxhint[d];
        min[d] = minhint[d];
    }

    /*
    printf("trysplit @ %g (%g %g %g %g %g %g) dim = %d, %td %td\n",
            node->split, 
            max[0], max[1], max[2], min[0],  min[1],  min[2],  
            dim, node->start, node->size);
    */
    ptrdiff_t p = node->start;

    double split;
    int dim = kd_build_find_longest_dim(Nd, min, max);

    /* recover from an imbalanced split */

    if(kd_build_is_poor_split(node, p)) {
        /* use the middle of the hinted box */
        split = (max[dim] + min[dim]) * 0.5;
        p = kd_build_pivot(node, dim, split);
    }

    if(kd_build_is_poor_split(node, p)) {
        /* Use the median of three samples as a split */
        double a = kd_input(tree, node->start, dim);
        double b = kd_input(tree, node->start + node->size / 2, dim);
        double c = kd_input(tree, node->start + node->size - 1, dim);
        /* trick from http://stackoverflow.com/a/19045659 */
        /* to prove assume a < b without loosing generality. */
        split = fmax(fmin(a, b), fmin(fmax(a, b), c));
        p = kd_build_pivot(node, dim, split);
    }

    if(kd_build_is_poor_split(node, p)) {
        /* shrink min max and retry, this is slow O(n),
         * thus we don't want to be here often */
        kd_build_update_min_max(node, min, max);

        /* use a new longest dim because min max are updated */
        dim = kd_build_find_longest_dim(Nd, min, max);

        split = (max[dim] + min[dim]) * 0.5;
        p = kd_build_pivot(node, dim, split);
    }

    if(kd_build_is_poor_split(node, p)) {
        /* if we are here, then there is no way to split the node. All
         * input are very close      */
        return next;
    }

    /* XXX: check the comment below one ineq shall be strict
     * after sliding we have data[0, .... p - 1] <= split
     * and data[q +1.... end] >= split */

    node->split = split;
    node->dim = dim;
    node->link[0] = kd_alloc(tree);
    node->link[0]->index = next++;
    node->link[0]->start = node->start;
    node->link[0]->size = p - node->start;
    node->link[0]->dim = -1;
    node->link[1] = kd_alloc(tree);
    node->link[1]->index = next++;
    node->link[1]->start = p;
    node->link[1]->size = node->size - (p - node->start);
    node->link[1]->dim = -1;
/*
    printf("will split %g (%td %td), (%td %td)\n", 
            *(double*)split, 
            node->link[0]->start, node->link[0]->size,
            node->link[1]->start, node->link[1]->size);
*/
    double mid[Nd];
    for(d = 0; d < Nd; d++) {
        mid[d] = max[d];
    }
    mid[node->dim] = node->split;
    next = kd_build_split(node->link[0], min, mid, next);

    for(d = 0; d < Nd; d++) {
        mid[d] = min[d];
    }
    mid[node->dim] = node->split;
    next = kd_build_split(node->link[1], mid, max, next);

    double * max1 = kd_node_max(node->link[1]);
    double * min1 = kd_node_min(node->link[1]);
    for(d = 0; d < Nd; d++) {
        max[d] = kd_node_max(node->link[0])[d];
        if(max[d] < max1[d]) max[d] = max1[d];
        min[d] = kd_node_min(node->link[0])[d];
        if(min[d] > min1[d]) min[d] = min1[d];
    }
    return next;
}

/* 
 * create a root KDNode based on input data specified in KDTree 
 * free it with kd_free
 * */
KDNode * 
kd_build(KDTree * tree) 
{
    ptrdiff_t i;
    int Nd = tree->input.dims[1];
    double min[Nd];
    double max[Nd];
    int d;
    if (tree->ind_size > 0) {
        for(d = 0; d < Nd; d++) {
            min[d] = kd_input(tree, 0, d);
            max[d] = kd_input(tree, 0, d);
        }
    }
    else {
        for(d = 0; d < Nd; d++) {
            min[d] = 0;
            max[d] = 0;
        }
    }
    for(i = 0; i < tree->ind_size; i++) {
        for(d = 0; d < Nd; d++) {
            double data = kd_input(tree, i, d);
            if(min[d] > data) { min[d] = data; }
            if(max[d] < data) { max[d] = data; }
        }
    }
    KDNode * root = kd_alloc(tree);
    root->start = 0;
    root->index = 0;
    root->dim = -1;
    root->size = tree->ind_size;
    tree->size = kd_build_split(root, min, max, 1);

    return root;
}

static void
kd_free0(KDTree * tree, size_t size, void * ptr) 
{
    if(tree->free == NULL) {
        free(ptr);
    } else {
        tree->free(tree->userdata, size, ptr);
    }
}

/**
 * free a tree from a node.
 * this is recursive
 * */
void 
kd_free(KDNode * node) 
{
    if(node->link[0]) kd_free(node->link[0]);
    if(node->link[1]) kd_free(node->link[1]);
    node->tree->size --;
    kd_free0(node->tree, 
            sizeof(KDNode) +
            sizeof(double) * 2 * node->tree->input.dims[1],
            node);
}

static double * 
kd_attr_init_r(KDAttr * attr, KDNode * node)
{
    double * rt = &attr->buffer[node->index * attr->input.dims[1]];
    ptrdiff_t d;
    for(d = 0; d < attr->input.dims[1]; d++) {
        rt[d] = 0;
    }
    if (node->dim < 0) {
        ptrdiff_t i;
        for(i = node->start; i < node->start + node->size; i ++) {
            for(d = 0; d < attr->input.dims[1]; d++) {
                rt[d] += kd_attr_get(attr, i, d);
            }
        }
        return rt;
    }
     
    double * left = kd_attr_init_r(attr, node->link[0]);
    double * right = kd_attr_init_r(attr, node->link[1]);
    for(d = 0; d < attr->input.dims[1]; d++) {
        rt[d] = left[d] + right[d];
    }
    return rt;
}

void 
kd_attr_init(KDAttr * attr, KDNode * root) 
{
    kd_attr_init_r(attr, root); 
}

