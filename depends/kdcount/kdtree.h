#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

typedef double (*kd_castfunc)(void * p);
typedef void (*kd_freefunc)(void* data, size_t size, void * ptr);
typedef void * (*kd_mallocfunc)(void* data, size_t size);

/* The following functions reutnr cull */
/* dist is a scalar output*/
typedef int (*kd_point_point_cullmetric)(void * userdata, int ndims, double * dx, double * dist);
/* distmin, distmax are scalar outputs*/
typedef int (*kd_node_node_cullmetric)(void * userdata, int ndims, double * min, double * max,
            double *distmin, double *distmax);

typedef struct KDArray {
    /* the buffer holding array elements required */
    char * buffer; 
    /* number of points. required*/
    ptrdiff_t dims[2];
    /* the byte offset of the axes.  required
     * the i-th position , d-th component is
     * at i * strides[0] + d * strides[1] */
    ptrdiff_t strides[2]; 

    /* if cast p1 to double and return it */
    double (* cast)(void * p1);
    /* the byte size of each scalar, required */
    ptrdiff_t elsize;

} KDArray;

typedef struct KDTree {

/* defining the input positions */
    KDArray input;
    /* a permutation array for indexing input, required and will be modified */
    ptrdiff_t * ind; 
    /* length of ind_size*/
    size_t ind_size;

/* the following defines how the tree is constructed */

    /* split thresh, required. 10 is good*/
    int thresh;
    /* the periodic boxsize per axis,
     * or NULL if there is no box  */
    double * boxsize;
    /* unused */
    double p;

/* memory allocation for KDNodes */
    /* allocate memory, NULL to use malloc() */
    kd_mallocfunc malloc;
    /* deallocate memory, size is passed in for a slab allocator,
     * NULL to use free() */
    kd_freefunc free;
    void * userdata;
    size_t size;
} KDTree;

typedef struct KDNode {
    KDTree * tree;
    ptrdiff_t index;
    struct KDNode * link[2];
    ptrdiff_t start;
    ptrdiff_t size;
    int dim; /* -1 for leaf */
    double split;
    char ext[];
} KDNode;

typedef struct KDAttr {
    KDTree * tree;
    KDArray input;

/* internal storage for cumulative weight on the nodes */
    double * buffer;
} KDAttr;

static inline double * kd_node_max(KDNode * node) {
    /* upper limit of the node */
    return (double*) (node->ext);
}
static inline double * kd_node_min(KDNode * node) {
    /* lower limit of the node */
    return kd_node_max(node) + node->tree->input.dims[1];
}
static inline double kd_array_get(KDArray * array, ptrdiff_t i, ptrdiff_t d) {
    char * ptr = & array->buffer[
                        i * array->strides[0] + 
                        d * array->strides[1]];
    if(array->cast) {
        return array->cast(ptr);
    } else {
        return * (double*) ptr;
    }
}

static inline double kd_input(KDTree * tree, ptrdiff_t i, ptrdiff_t d) {
    i = tree->ind[i];
    return kd_array_get(&tree->input, i, d);
}
static inline void kd_swap(KDTree * tree, ptrdiff_t i, ptrdiff_t j) {
    ptrdiff_t t = tree->ind[i];
    tree->ind[i] = tree->ind[j];
    tree->ind[j] = t;
}

static inline double * kd_attr_get_node(KDAttr * attr, KDNode * node) {
    return &attr->buffer[node->index * attr->input.dims[1]];
}

static inline double kd_attr_get(KDAttr * attr, ptrdiff_t i, ptrdiff_t d) {
    return kd_array_get(&attr->input, attr->tree->ind[i], d);
}

static inline void
kd_realminmax(KDTree * tree, double min, double max, double * realmin, double * realmax, int d) 
{
    if(tree->boxsize) {
        double full = tree->boxsize[d];
        double half = full * 0.5;
        /* periodic */
        /* /\/\ */
        if(max <= 0 || min >= 0) {
            /* do not pass through 0 */
            min = fabs(min);
            max = fabs(max);
            if(min > max) {
                double t = min;
                min = max;
                max = t;
            }
            if(max < half) {
                /* all below half*/
                *realmin = min;
                *realmax = max;
            } else if(min > half) {
                /* all above half */
                *realmax = full - min;
                *realmin = full - max;
            } else {
                /* min below, max above */
                *realmax = half;
                *realmin = fmin(min, full - max);
            }
        } else {
            /* pass though 0 */
            min = -min;
            if(min > max) max = min;
            if(max > half) max = half;
            *realmax = max;
            *realmin = 0;
        }
    } else {
        /* simple */
        /* \/     */
        if(max <= 0 || min >= 0) {
            /* do not pass though 0 */
            min = fabs(min);
            max = fabs(max);
            if(min < max) {
                *realmin = min;
                *realmax = max;
            } else {
                *realmin = max;
                *realmax = min;
            }
        } else {
            min = fabs(min);
            max = fabs(max);
            *realmax = fmax(max, min);
            *realmin = 0;
        }
    }

}
static inline double
kd_realdiff(KDTree * tree, double dx, int d) 
{
    if (dx < 0) dx = - dx;
    if (tree->boxsize) {
        double half = 0.5 * tree->boxsize[d];
        if (dx > half) dx = tree->boxsize[d] - dx;
    }
    return dx;
}
static inline void kd_collect(KDNode * node, KDArray * input, double * ptr) {

    /* collect permuted elements into a double array, 
     * so that they can be paired quickly (cache locality!)*/
    ptrdiff_t * ind = node->tree->ind;
    int Nd = input->dims[1];
    char * base = input->buffer;
    ptrdiff_t j;
    for (j = node->start; j < node->start + node->size; j++) {
        int d;
        char * item = base + ind[j] * input->strides[0];
        if(input->cast) {
            for(d = 0; d < Nd; d++) {
                *ptr = input->cast(item);
                ptr++;
                item += input->strides[1];
            }
        } else {
            for(d = 0; d < Nd; d++) {
                memcpy(ptr, item, input->elsize);
                ptr++;
                item += input->strides[1];
            }
        }
    }
}

KDNode *
kd_build(KDTree * tree);

void
kd_free(KDNode * node);

void 
kd_attr_init(KDAttr * attr, KDNode * root);

/* x and y points to the double position */
typedef struct KDEnumPair {
    double r;
    ptrdiff_t i;
    ptrdiff_t j;
} KDEnumPair;

typedef struct KDEnumNodePair {
    double distmax2;
    double distmin2;
    KDNode * nodes[2];
} KDEnumNodePair;

/* returns whether the edge enumeration between these two nodes shall be terminated */
typedef int (*kd_enum_visit_edge)(void * userdata, KDEnumPair * pair);
/* returns whether the edge enumeration shall be terminated */
typedef int (*kd_enum_check_nodes)(void * userdata, KDEnumNodePair * pair);
/* returns whether the child nodes, if exist, shall be opened. */
typedef int (*kd_enum_visit_node)(void * userdata, KDNode * node);

int
kd_enum_check(KDNode * nodes[2], double maxr2, int skip_symmetric, kd_enum_visit_edge visit_edge, void * userdata);

int
kd_enum(KDNode * nodes[2], double maxr,
        kd_enum_visit_edge visit_edge,
        kd_enum_check_nodes check_nodes,
        void * userdata);

int
kd_enum_full(KDNode * nodes[2], double maxr,
        kd_enum_visit_edge visit_edge,
        kd_enum_check_nodes check_nodes,
        kd_enum_visit_node visit_node,
        double opening_factor,
        int skip_symmetric,
        void * userdata);

int
kd_fof(KDNode * tree, double linking_length, ptrdiff_t * head);
int
kd_fof_linkedlist(KDNode * tree, double linking_length, ptrdiff_t * head);
int
kd_fof_allpairs(KDNode * tree, double linking_length, ptrdiff_t * head);
int 
kd_fof_unsafe(KDNode * tree, double linking_length, ptrdiff_t * head);
int 
kd_fof_heuristics(KDNode * tree, double linking_length, ptrdiff_t * head);

void
kd_count(KDNode * nodes[2], KDAttr * attrs[2],
        double * edges, uint64_t * count, double * weight,
        int nedges,
        kd_point_point_cullmetric ppcull,
        kd_node_node_cullmetric nncull,
        void * userdata,
        uint64_t * brute_force,
        uint64_t * node_node
        );

void
kd_integrate(KDNode * node, KDAttr * attr,
        uint64_t * count, double * weight,
        double * min, double * max,
        uint64_t * brute_force,
        uint64_t * node_node
);

/* function called to compute the force between two center of masses;
 * f is force / (m1 * m2). Thus it makes sense only for 'additive' force
 * like gravity or pair counting.
 * */

typedef void (*kd_force_func)(double r, double * dx, double * f, int ndims, void * userdata);

/*
 * compute the force between pos and all particles in a  tree node.
 * mass and xmass are mass and mass weighted position of particles in the tree ndoe.
 *
 * r_cut is the cut-off scale.
 * eta is the opening criteria for the tree force. set eta to 0 to use exact pair-wise force.
 * force is an output array of size ndim.
 *
 * Currently the dimension of force and the dimension of position must be the same.
 * */
void
kd_force(double * pos, KDNode * node, KDAttr * mass, KDAttr * xmass,
        double r_cut, double eta, double * force,
        kd_force_func func, void * userdata);
