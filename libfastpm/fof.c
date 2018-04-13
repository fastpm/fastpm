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
    int ThisTask;
    int NTask;
};

/* creating a kdtree struct
 * for store with np particles starting from start.
 * */
static 
KDNode *
_create_kdtree (KDTree * tree, FastPMStore * store, ptrdiff_t start, size_t np, double boxsize[])
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


    tree->thresh = 8;

    tree->boxsize = boxsize;

    KDNode * root = kd_build(tree);
    fastpm_info("Creating KDTree with %td nodes for %td particles\n", tree->size, np);
    return root;
}
static void
_update_local_minid(FastPMFOFFinder * finder, 
    struct FastPMFOFData * fof, uint64_t * id, ptrdiff_t * head, size_t np, size_t nghosts)
{
    /* update min id like a hash table */
    ptrdiff_t i;
    /* initialize minid, used as a global tag of groups as we merge */
    for(i = 0; i < np; i ++) {
        fof[i].minid = id[i];
        fof[i].task = finder->priv->ThisTask;
    }

    for(i = np; i < np + nghosts; i ++) {
        fof[i].minid = id[i];
        fof[i].task = -1;
    }

    /* all of the unique heads are updated. first */
    for(i = 0; i < np + nghosts; i ++) {
        if(id[i] < fof[head[i]].minid) {
            fof[head[i]].minid = id[i];
            /* if minid is a ghost the hosting rank is undetermined. */
            fof[head[i]].task = (i < np) ? finder->priv->ThisTask : -1;
        }
    }
}

static void
_sync_local_minid(struct FastPMFOFData * fof, ptrdiff_t * head, size_t np, size_t np_ghosts)
{
    ptrdiff_t i;
    /* then sync the non head items */
    for(i = 0; i < np + np_ghosts; i++) {
        fof[i].minid = fof[head[i]].minid;
        fof[i].task = fof[head[i]].task;
    }
}
void
fastpm_fof_init(FastPMFOFFinder * finder, FastPMStore * store, PM * pm)
{
    finder->priv = malloc(sizeof(FastPMFOFFinderPrivate));
    finder->p = store;
    finder->pm = pm;
    MPI_Comm_rank(pm_comm(pm), &finder->priv->ThisTask);
    MPI_Comm_size(pm_comm(pm), &finder->priv->NTask);
}

/*
static int
_visit_edge(void * userdata, KDEnumPair * pair)
{
    return 0;
}
*/
static void
_reduce_minid(PMGhostData * pgd, enum FastPMPackFields attributes,
        ptrdiff_t i, void * buffer, void * userdata)
{
    void ** data = (void**) userdata;
    ptrdiff_t * head = data[0];
    ptrdiff_t * merged = data[1];
    int * iter = data[2];

    int ThisTask;
    MPI_Comm_rank(pm_comm(pgd->pm), &ThisTask);

    FastPMStore * p = pgd->p;
    struct FastPMFOFData * remote = buffer;

    int merge = 0;

    {
        ptrdiff_t i;
        for(i = 0; i < p->np; i ++) {
            if(p->fof[head[i]].minid > p->id[i]) {
                fastpm_raise(-1, "raise = %d %d merge %d id = %03d minid %03td -> %03td task %03d -> %03d merge = %d nmerged = %d\n", i, head[i], ThisTask, p->id[i], remote->minid, p->fof[head[i]].minid, remote->task, p->fof[head[i]].task, merge, *merged);

            }
        }
    }

    if(remote->task >= 0) {
        /* minid is on that processor */
        if(p->fof[head[i]].task == -1
        /* but not on this processor */
        ) {
             merge = 1;
        }
        /* or min id reduces */
        if(remote->minid < p->fof[head[i]].minid) {
            merge = 1;
        }
    }

    if(p->id[i] == 198) {
        fastpm_ilog(INFO, "iter = %d merge %d id = %03d minid %03td -> %03td task %03d -> %03d merge = %d nmerged = %d\n", *iter, ThisTask, p->id[i], remote->minid, p->fof[head[i]].minid, remote->task, p->fof[head[i]].task, merge, *merged);
    }

    if(merge) {
        p->fof[head[i]].minid = remote->minid;
        p->fof[head[i]].task = remote->task;

        (*merged) ++;
    }
}


static int
FastPMTargetFOF(FastPMStore * store, ptrdiff_t i, void ** userdata)
{
    struct FastPMFOFData * fof = userdata[0];
    ptrdiff_t * head = userdata[1];
    return fof[head[i]].task;
}

static void
_send(PMGhostData * pgd, struct FastPMFOFData * fof, ptrdiff_t *head)
{
    FastPMStore * p = pgd->p;

    struct FastPMFOFData * fof_ghosts = malloc(sizeof(fof[0]) * pgd->nghosts);

    memcpy(fof_ghosts, &fof[p->np], sizeof(fof[0]) * pgd->nghosts);

    pm_ghosts_send(pgd);

    ptrdiff_t i;
    for(i = p->np; i < p->np + pgd->nghosts; i ++) {
        int restore = 0;
        if(fof_ghosts[i - p->np].minid < fof[head[i]].minid ) {
            restore = 1;
        }

        if(fof[i].task == -1){
            restore = 1;
        }

        if (restore) {
            fof[head[i]].minid = fof_ghosts[i - p->np].minid;
            fof[head[i]].task = fof_ghosts[i - p->np].task;
        }
    }

    free(fof_ghosts);
}

static void
fastpm_fof_decompose(FastPMFOFFinder * finder, FastPMStore * p, PM * pm)
{
    struct FastPMFOFData * fof = malloc(sizeof(fof[0]) * p->np_upper);
    ptrdiff_t * head = malloc(sizeof(head[0]) * p->np_upper);
    uint64_t * id = p->id;
    ptrdiff_t i;
    
    p->fof = fof;

    /* initial decompose -- reduce number of ghosts */
    fastpm_store_decompose(p,
                (fastpm_store_target_func) FastPMTargetPM,
                pm, pm_comm(pm));
    
    /* create ghosts mesh size is usually > ll so we are OK here. */
    double below[3], above[3];

    int d;
    for(d = 0; d < 3; d ++) {
        /* bigger padding reduces number of iterations */
        below[d] = -finder->linkinglength;
        above[d] = finder->linkinglength;
    }

    /* merge */
    PMGhostData * pgd = pm_ghosts_create_full(pm, p,
            PACK_POS | PACK_ID | PACK_FOF, NULL,
            below, above
        );


    KDTree tree;

    KDNode * root = _create_kdtree(&tree, p, 0, p->np + pgd->nghosts, pm_boxsize(pm));

    /* local find */

    kd_fof(root, finder->linkinglength, head);

    _update_local_minid(finder, fof, id, head, p->np, pgd->nghosts);

    {
        int iter = 0;
        while(1) {
            size_t nmerged = 0;
            void * userdata[3] = {head, &nmerged, &iter};

            _sync_local_minid(fof, head, p->np, pgd->nghosts);

            _send(pgd, fof, head);

            for(i = 0; i < p->np + pgd->nghosts; i ++) {
                if(p->id[i] == 198)
                fastpm_ilog(INFO, "iter%d %d id = %03td head=%d, head[%03d/%03d] : id = %03td task = %03d\n", iter, finder->priv->ThisTask, p->id[i], i, p->id[head[i]], p->np, fof[i].minid, fof[i].task);
            }

            _sync_local_minid(fof, head, p->np, pgd->nghosts);

            pm_ghosts_reduce_any(pgd, PACK_FOF, _reduce_minid, userdata);

            for(i = 0; i < p->np + pgd->nghosts; i ++) {
                if(p->id[i] == 198)
                fastpm_ilog(INFO, "iter%d %d id = %03td head=%d, head[%03d/%03d] : id = %03td task = %03d\n", iter, finder->priv->ThisTask, p->id[i], i, p->id[head[i]], p->np, fof[i].minid, fof[i].task);
            }

            MPI_Allreduce(MPI_IN_PLACE, &nmerged, 1, MPI_LONG, MPI_SUM, pm_comm(pm));

            MPI_Barrier(pm_comm(pm));

            fastpm_info("FOF reduction iteration %d : merged %td crosslinks\n", iter, nmerged);

            if(nmerged == 0) break;
            iter++;
        }
    }

    p->fof = NULL;

    pm_ghosts_free(pgd);

    size_t n_nonlocal = 0;
    for(i = 0; i < p->np; i ++) {
        if(fof[i].task == -1) {
            printf("check %d id = %03td head = %03td head[%03d/%03d] : id = %03td task = %03d\n", finder->priv->ThisTask, p->id[i], p->id[head[i]], i, p->np, fof[i].minid, fof[head[i]].task);
        }
        if(fof[i].task != finder->priv->ThisTask) {
            n_nonlocal ++;
        }
        //printf("%d: %ld fof = {%ld, %d} \n", finder->priv->ThisTask, p->id[i], p->fof[i].minid, p->fof[i].task);
    }
    fastpm_info("%td particles are linked to a remote group\n", n_nonlocal);

    /* real decompose, move particles of the same halo to the same rank */
    {
        void * userdata[2] = {fof, head};
        /* attributes */
        fastpm_store_decompose(p,
                    (fastpm_store_target_func) FastPMTargetFOF,
                    userdata, pm_comm(pm));
    }

    free(head);
    free(fof);
    kd_free(root);

    free(tree.ind);
}

static size_t
_assign_halo_attr(ptrdiff_t * head, ptrdiff_t * offset, size_t np, int nmin)
{
    ptrdiff_t i;
    for(i = 0; i < np; i ++) {
        offset[i] = 0;
    }

    /* set offset to number of particles in the halo */
    for(i = 0; i < np; i ++) {
        offset[head[i]] ++;
    }

    size_t it = 0;
    
    /* assign attr index */
    for(i = 0; i < np; i ++) {
        if(offset[i] >= nmin) {
            offset[i] = it;
            it ++;
        } else {
            offset[i] = -1;
        }
    }
    return it;
}

void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos)
{
    fastpm_fof_decompose(finder, finder->p, finder->pm);

    KDTree tree;

    /* redo fof on the new decomposition -- no halo cross two ranks */
    KDNode * root = _create_kdtree(&tree, finder->p, 0, finder->p->np, pm_boxsize(finder->pm));

    ptrdiff_t * head = malloc(sizeof(head[0]) * finder->p->np);
    ptrdiff_t * offset = malloc(sizeof(offset[0]) * finder->p->np);

    kd_fof(root, finder->linkinglength, head);

    /* now do the attributes */

    /* */
    {
        size_t nhalos = _assign_halo_attr(head, offset, finder->p->np, finder->nmin);
    /*
        ptrdiff_t i;
        for(i = 0; i < finder->p->np; i++) {
            if(head[i] < 0) abort();
            if(head[i] >= finder->p->np) abort();
        }
*/
        fastpm_store_init(halos, nhalos, (finder->p->attributes | PACK_LENGTH) & ~PACK_ID
                , FASTPM_MEMORY_HEAP);
        halos->np = nhalos;

        MPI_Allreduce(MPI_IN_PLACE, &nhalos, 1, MPI_LONG, MPI_SUM, pm_comm(finder->pm));

        fastpm_info("Found %td halos >= %d particles\n", nhalos, finder->nmin);
    }

    ptrdiff_t i;
    for(i = 0; i < halos->np; i++) {
        int d;
        halos->length[i] = 0;

        if(halos->aemit)
            halos->aemit[i] = 0;

        for(d = 0; d < 3; d++) {
            if(halos->x)
                halos->x[i][d] = 0;
            if(halos->v)
                halos->v[i][d] = 0;
            if(halos->dx1)
                halos->dx1[i][d] = 0;
            if(halos->dx2)
                halos->dx2[i][d] = 0;
            if(halos->q)
                halos->q[i][d] = 0;
        }
    }
    for(i = 0; i < finder->p->np; i++) {
        ptrdiff_t hid = offset[head[i]];
        if(hid < 0) continue;

        if(hid >= halos->np) {
            abort();
        }
        halos->length[hid] += 1;

        if(halos->aemit)
            halos->aemit[hid] += finder->p->aemit[i];
        int d;
        for(d = 0; d < 3; d++) {
            if(halos->x)
                halos->x[hid][d] += finder->p->x[i][d];
            if(halos->v)
                halos->v[hid][d] += finder->p->v[i][d];
            if(halos->dx1)
                halos->dx1[hid][d] += finder->p->dx1[i][d];
            if(halos->dx2)
                halos->dx2[hid][d] += finder->p->dx2[i][d];
            if(halos->q)
                halos->q[hid][d] += finder->p->q[i][d];
        }
    }

    for(i = 0; i < halos->np; i++) {
        int d;
        double n = halos->length[i];

        for(d = 0; d < 3; d++) {
            if(halos->x)
                halos->x[i][d] /= n;
            if(halos->v)
                halos->v[i][d] /= n;
            if(halos->dx1)
                halos->dx1[i][d] /= n;
            if(halos->dx2)
                halos->dx2[i][d] /= n;
            if(halos->q)
                halos->q[i][d] /= n;
        }
        if(halos->aemit)
            halos->aemit[i] /= n;

    }

    free(head);
    kd_free(root);
    free(tree.ind);
}

void
fastpm_fof_destroy(FastPMFOFFinder * finder)
{
    free(finder->priv);
}


