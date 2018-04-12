#include <string.h>
#include <alloca.h>
#include <math.h>
#include <mpi.h>

#include <kdcount/kdtree.h>

#include <fastpm/libfastpm.h>
#include <gsl/gsl_rng.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/store.h>

#include <fastpm/fof.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "pmghosts.h"
#include "vpm.h"

struct FastPMFOFFinderPrivate {
    int foo;
};

/* creating a kdtree struct
 * for store with np particles starting from start.
 * */
static 
KDNode *
_create_kdtree (KDTree * tree, FastPMStore * store, ptrdiff_t start, size_t np)
{
    ptrdiff_t * ind = malloc(sizeof(ind[0]) * np);
    ptrdiff_t i;

    for(i = 0; i < np; i ++) {
        ind[i] = start + i;
    }

    tree->input.buffer = (char*) store->x;
    tree->input.dims[0] = np;
    tree->input.dims[1] = 3;
    tree->input.strides[0] = sizeof(store->x[0]);
    tree->input.strides[1] = sizeof(store->x[0][0]);
    tree->input.elsize = sizeof(store->x[0][0]);
    tree->input.cast = NULL;

    tree->ind = ind;
    tree->ind_size = np;
    tree->malloc = NULL;
    tree->free = NULL;


    tree->thresh = 10;

    tree->boxsize = NULL;

    KDNode * root = kd_build(tree);
    fastpm_info("Creating KDTree with %td nodes for %td particles\n", tree->size, np);
    return root;
}
static void
_update_local_minid(uint64_t * minid, uint64_t * id, ptrdiff_t * head, size_t np)
{
    /* update min id; */

    ptrdiff_t i;
    /* initialize minid, used as a global tag of groups as we merge */
    for(i = 0; i < np; i ++) {
        minid[i] = id[i];
    }

    /* all of the unique heads are updated. first */
    for(i = 0; i < np; i ++) {
        if(id[i] < minid[head[i]]) {
            minid[head[i]] = id[i];
        }
    }
    /* then copy the minid of my head */
    for(i = 0; i < np; i++) {
        minid[i] = minid[head[i]];
    }
}

void
fastpm_fof_init(FastPMFOFFinder * finder, FastPMStore * store, PM * pm)
{
    finder->priv = malloc(sizeof(FastPMFOFFinderPrivate));
    finder->p = store;
    finder->pm = pm;
}
static int
_visit_edge(void * userdata, KDEnumPair * pair)
{
    return 0;
}

static void
myreduction(PMGhostData * pgd, enum FastPMPackFields attributes,
        ptrdiff_t index, void * buffer, void * userdata)
{

}

void
fastpm_fof_execute(FastPMFOFFinder * finder)
{
    uint64_t * minid = malloc(sizeof(uint64_t ) * finder->p->np_upper);
    ptrdiff_t * head = malloc(sizeof(ptrdiff_t) * finder->p->np_upper);
    uint64_t * id = finder->p->id;
    //ptrdiff_t i;
    
    FastPMStore * p = finder->p;
    p->fof_minid = minid;

    /* decompose */
    fastpm_store_decompose(p,
                (fastpm_store_target_func) FastPMTargetPM,
                finder->pm, pm_comm(finder->pm));
    
    /* create ghosts mesh size is usually > ll so we are OK here. */
    PMGhostData * pgd = pm_ghosts_create(finder->pm, p,
            PACK_POS | PACK_ID | PACK_FOF_MINID, NULL);

    KDTree tree, tree_ghosts;

    KDNode * root = _create_kdtree(&tree, p, 0, p->np + pgd->nghosts);
    KDNode * root_ghosts = _create_kdtree(&tree_ghosts, p, p->np, pgd->nghosts);

    /* local find */

    kd_fof(root, finder->linkinglength, head);
    
    _update_local_minid(minid, id, head, p->np + pgd->nghosts);

    pm_ghosts_reduce_any(pgd, PACK_FOF_MINID, myreduction, NULL);
    /* merge */
    KDNode * nodes[2];
    nodes[0] = root_ghosts;
    nodes[1] = root;

    kd_enum(nodes, finder->linkinglength, _visit_edge, NULL, NULL);

    /* attributes */
    
    pm_ghosts_free(pgd);

    free(minid);
    free(head);
    p->fof_minid = NULL;
    kd_free(root);
    kd_free(root_ghosts);

    free(tree.ind);
    free(tree_ghosts.ind);
}

void
fastpm_fof_destroy(FastPMFOFFinder * finder)
{
    free(finder->priv);
}


