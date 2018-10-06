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

#include <fastpm/io.h>
#include <fastpm/string.h>
#include <bigfile-mpi.h>

struct FastPMFOFFinderPrivate {
    int ThisTask;
    int NTask;
    double * boxsize;
};

/* creating a kdtree struct
 * for store with np particles starting from start.
 * */
struct KDTreeNodeBuffer {
    void * mem;
    void * base;
    char * ptr;
    char * end;
    struct KDTreeNodeBuffer * prev;
};

static void
fastpm_fof_compute_local_halo_attrs (FastPMFOFFinder * finder, FastPMStore * halos, struct FastPMFOFData * fofsave, ptrdiff_t * head, size_t np_and_ghosts);

static void
fastpm_fof_reduce_sorted_local_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos);

static void
fastpm_fof_add_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos, int hid, FastPMStore * h1, int i);

static void
fastpm_fof_create_local_halos(FastPMFOFFinder * finder, FastPMStore * halos, size_t nhalos);

static void
fastpm_fof_convert_particle_to_halo(FastPMStore * p, ptrdiff_t i, FastPMStore * halos);

static void *
_kdtree_buffered_malloc(void * userdata, size_t size)
{
    struct KDTreeNodeBuffer ** pbuffer = (struct KDTreeNodeBuffer**) userdata;

    struct KDTreeNodeBuffer * buffer = *pbuffer;

    if(buffer->base == NULL || buffer->ptr + size >= buffer->end) {
        struct KDTreeNodeBuffer * newbuffer = malloc(sizeof(newbuffer[0]));
        newbuffer->mem = buffer->mem;
        size_t newsize = 1024 * 1024 * 32; /* 32 MB for each block */
        if(newsize < size) {
            newsize = size;
        }
        newbuffer->base = fastpm_memory_alloc(buffer->mem, "KDTreeBase", newsize, FASTPM_MEMORY_STACK);
        newbuffer->ptr = newbuffer->base;
        newbuffer->end = newbuffer->base + newsize;
        newbuffer->prev = buffer;

        *pbuffer = newbuffer;
        buffer = newbuffer;
    }

    void * r = buffer->ptr;
    buffer->ptr += size;
    return r;
}

static void
_kdtree_buffered_free(void * userdata, size_t size, void * ptr)
{
    /* do nothing; */
}

static 
KDNode *
_create_kdtree (KDTree * tree, int thresh,
    FastPMStore ** stores, int nstore,
    double boxsize[])
{
    /* if boxsize is NULL the tree will be non-periodic. */
    /* the allocator; started empty. */
    struct KDTreeNodeBuffer ** pbuffer = malloc(sizeof(void*));
    struct KDTreeNodeBuffer * headbuffer = malloc(sizeof(headbuffer[0]));
    *pbuffer = headbuffer;

    headbuffer->mem = stores[0]->mem;
    headbuffer->base = NULL;
    headbuffer->prev = NULL;

    int s;
    ptrdiff_t i;

    tree->userdata = pbuffer;

    tree->input.dims[0] = 0;

    for(s = 0; s < nstore; s ++) {
        tree->input.dims[0] += stores[s]->np;
    }

    tree->input.dims[1] = 3;

    if(tree->input.dims[0] < stores[0]->np_upper) {
        /* if the first store is big enough, use it for the tree */
        tree->input.buffer = (void*) &stores[0]->x[0][0];
    } else {
        /* otherwise, allocate a big buffer and make a copy */
        tree->input.buffer = _kdtree_buffered_malloc(pbuffer,
                    tree->input.dims[0] * sizeof(stores[0]->x[0]));
        memcpy(tree->input.buffer, &stores[0]->x[0][0], stores[0]->np * sizeof(stores[0]->x[0]));
    }

    i = stores[0]->np;

    /* copy the other positions to the base pointer. */
    for(s = 1; s < nstore; s ++) {
        memcpy(((char*) tree->input.buffer) + i * sizeof(stores[0]->x[0]),
                &stores[s]->x[0][0],
                stores[s]->np * sizeof(stores[0]->x[0]));

        i = i + stores[s]->np;
    }

    tree->input.strides[0] = sizeof(stores[0]->x[0]);
    tree->input.strides[1] = sizeof(stores[0]->x[0][0]);
    tree->input.elsize = sizeof(stores[0]->x[0][0]);
    tree->input.cast = NULL;

    tree->ind = _kdtree_buffered_malloc(pbuffer, tree->input.dims[0] * sizeof(tree->ind[0]));
    for(i = 0; i < tree->input.dims[0]; i ++) {
        tree->ind[i] = i;
    }
    tree->ind_size = tree->input.dims[0];

    tree->malloc = _kdtree_buffered_malloc;
    tree->free = _kdtree_buffered_free;

    tree->thresh = thresh;

    tree->boxsize = boxsize;

    KDNode * root = kd_build(tree);
    fastpm_info("Creating KDTree with %td nodes for %td particles\n", tree->size, tree->ind_size);
    return root;
}

void 
_free_kdtree (KDTree * tree, KDNode * root)
{
    kd_free(root);
    struct KDTreeNodeBuffer * buffer, *q, **pbuffer = tree->userdata;

    for(buffer = *pbuffer; buffer; buffer = q) {
        if(buffer->base)
            fastpm_memory_free(buffer->mem, buffer->base);
        q = buffer->prev;
        free(buffer);
    }
    free(pbuffer);
}

void
fastpm_fof_init(FastPMFOFFinder * finder, FastPMStore * store, PM * pm)
{
    finder->priv = malloc(sizeof(FastPMFOFFinderPrivate));
    finder->p = store;
    finder->pm = pm;

    finder->event_handlers = NULL;

    if (finder->periodic)
        finder->priv->boxsize = pm_boxsize(pm);
    else
        finder->priv->boxsize = NULL;

    MPI_Comm_rank(pm_comm(pm), &finder->priv->ThisTask);
    MPI_Comm_size(pm_comm(pm), &finder->priv->NTask);
}

static void
_fof_local_find(FastPMFOFFinder * finder,
            FastPMStore * p,
            PMGhostData * pgd,
            ptrdiff_t * head, double linkinglength)
{
    /* local find of p and the the ghosts */
    KDTree tree;

    FastPMStore * stores[2] = {p, pgd->p};

    KDNode * root = _create_kdtree(&tree, finder->kdtree_thresh, stores, 2, finder->priv->boxsize);

    kd_fof(root, linkinglength, head);

    _free_kdtree(&tree, root);
}

static int
_merge(struct FastPMFOFData * remote, struct FastPMFOFData * fofhead)
{
    int merge = 0;
    if(remote->minid < fofhead->minid) {
        merge = 1;
    }
    if(remote->minid == fofhead->minid &&
       remote->task != fofhead->task &&
       remote->task >= 0) {
        merge = 1;
    }

    if(merge) {
        fofhead->minid = remote->minid;
        fofhead->task = remote->task;

    }
    return merge;
}

static int
_merge2(struct FastPMFOFData * remote, struct FastPMFOFData * fofhead)
{
    int merge = 0;
    if(remote->minid < fofhead->minid) {
        merge = 1;
    }
    if(remote->minid == fofhead->minid &&
       remote->task != fofhead->task &&
       remote->task >= 0) {
        merge = 1;
    }

    if(merge) {
/*
        if(fofhead->minid != remote->minid)
            printf("merging minid %ld setting to %ld\n", fofhead->minid, remote->minid);
*/
        fofhead->minid = remote->minid;
        fofhead->task = remote->task;
    }
    return merge;
}

static void
_fof_global_merge(
    FastPMFOFFinder * finder,
    PMGhostData * pgd,
    struct FastPMFOFData * fofsave,
    ptrdiff_t * head
)
{
    ptrdiff_t i;
    FastPMStore * p = finder->p;
    PM * pm = finder->pm;

    struct FastPMFOFData * fofcomm = fastpm_memory_alloc(p->mem, "FOFComm",
                    sizeof(fofcomm[0]) * (p->np + pgd->p->np), FASTPM_MEMORY_STACK);

    memset(fofcomm, 0, sizeof(fofcomm[0]) * (p->np + pgd->p->np));

    /* initialize minid, used as a global tag of groups as we merge */
    for(i = 0; i < p->np; i ++) {
        fofsave[i].minid = p->id[i];
        fofsave[i].task = finder->priv->ThisTask;
    }
    for(i = 0; i < pgd->p->np; i ++) {
        fofsave[i + p->np].minid = pgd->p->id[i];
        fofsave[i + p->np].task = -1;
    }

    /* reduce the minid of the head items according to the local connection. */

    for(i = 0; i < p->np + pgd->p->np; i ++) {
        _merge(&fofsave[i], &fofsave[head[i]]);
    }

    /* at this point all items with head[i] = i have local minid and task */

    int iter = 0;

    while(1) {

        /* prepare the communication buffer, every particle has
         * the minid and task of the current head. such that
         * they will connect to the other ranks correctly */
        for(i = 0; i < p->np + pgd->p->np; i ++) {
            fofcomm[i].minid = fofsave[head[i]].minid;
            fofcomm[i].task = fofsave[head[i]].task;
        }

        /* at this point all items on fofcomm have local minid and task, ready to send */

        p->fof = fofcomm;
        pm_ghosts_send(pgd, PACK_FOF);

        p->fof = NULL;

        /* flatten the fof data received from the ghosts */
        /* FIXME: use reduce? */ 
        for(i = p->np; i < p->np + pgd->p->np; i ++) {
            fofcomm[i] = pgd->p->fof[i - p->np];
        }
        size_t nmerged = 0;
        for(i = p->np; i < p->np + pgd->p->np; i ++) {
            int m = _merge2(&fofcomm[i], &fofsave[head[i]]);
            nmerged += m;
        }

        /* at this point all items on fofsave with head[i] = i have present minid and task */

        MPI_Allreduce(MPI_IN_PLACE, &nmerged, 1, MPI_LONG, MPI_SUM, pm_comm(pm));

        MPI_Barrier(pm_comm(pm));

        fastpm_info("FOF reduction iteration %d : merged %td crosslinks\n", iter, nmerged);

        if(nmerged == 0) break;

        for(i = 0; i < p->np + pgd->p->np; i ++) {
            if(fofsave[i].minid < fofsave[head[i]].minid) {
                fastpm_raise(-1, "fofsave invariance is broken i = %td np = %td\n", i, p->np);
            }
        }

        iter++;
    }

    for(i = 0; i < p->np; i ++) {
        if(fofcomm[i].task == -1) {
            fastpm_raise(-1, "undeterimined particle %d id = %03td head = %03td head[%03d/%03d] : id = %03td task = %03d\n",
                finder->priv->ThisTask, p->id[i], p->id[head[i]], i, p->np, fofcomm[i].minid, fofcomm[i].task);
        }
        fofsave[i].minid = fofcomm[i].minid;
        fofsave[i].task = fofcomm[i].task;
    }
    for(i = p->np; i < p->np + pgd->p->np; i ++) {
        fofsave[i].task = -1;
    }

    fastpm_memory_free(p->mem, fofcomm);
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
    
    ptrdiff_t max = 0;
    /* assign attr index */
    for(i = 0; i < np; i ++) {
        if(offset[i] > max) {
            max = offset[i];
        }
        if(offset[i] >= nmin) {
            offset[i] = it;
            it ++;
        } else {
            offset[i] = -1;
        }
    }

    size_t nhalos = it;

    /* update head [i] to offset[head[i]], which stores the index in the halo store for this particle. */
    for(i = 0; i < np; i ++) {
        head[i] = offset[head[i]];
        if(head[i] != -1 && head[i] > nhalos) {
            fastpm_raise(-1, "head[i] (%td) > nhalos (%td) This shall not happend\n", head[i], nhalos);
        }
    }

    return it;
}
static double
periodic_add(double x, double wx, double y, double wy, double L)
{
    if(wx > 0) {
        while(y - x > L / 2) y -= L;
        while(y - x < -L / 2) y += L;

        return wx * x + wy * y;
    } else {
        return wy * y;
    }
}

/* first every undecided halo; then the real halo with particles */
static int
FastPMLocalSortByMinID(const int i1,
                    const int i2,
                    FastPMStore * p)
{
    int d1 = (p->fof[i1].task == -1);
    int d2 = (p->fof[i2].task == -1);

    if(d1 != d2) return d2 - d1;

    int v1 = (p->fof[i1].minid < p->fof[i2].minid);
    int v2 = (p->fof[i1].minid > p->fof[i2].minid);

    return v2 - v1;
}

static int
FastPMTargetMinID(FastPMStore * store, ptrdiff_t i, void * userdata)
{
    FastPMFOFFinder * finder = userdata;

    if(store->fof[i].task != -1) {
        return store->fof[i].minid % finder->priv->NTask;
    } else {
        return finder->priv->ThisTask;
    }
}

static int
FastPMTargetTask(FastPMStore * store, ptrdiff_t i, void * userdata)
{
    FastPMFOFFinder * finder = userdata;

    if(store->fof[i].task != -1)
        return store->fof[i].task;
    else
        return finder->priv->ThisTask;
}


void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos)
{
    /* initial decompose -- reduce number of ghosts */
    FastPMStore * p = finder->p;
    PM * pm = finder->pm;

    /* only do wrapping for periodic data */
    if(finder->priv->boxsize)
        fastpm_store_wrap(p, finder->priv->boxsize);

    double npmax, npmin, npstd, npmean;

    MPIU_stats(pm_comm(pm), p->np, "<->s", &npmin, &npmean, &npmax, &npstd);

    fastpm_info("load balance after decompose : min = %g max = %g mean = %g std = %g\n",
        npmin, npmax, npmean, npstd
        );

    /* still route particles to the pm pencils as if they are periodic. */
    if(0 != fastpm_store_decompose(p,
                (fastpm_store_target_func) FastPMTargetPM,
                pm, pm_comm(pm))
    ) {
        fastpm_raise(-1, "out of storage space decomposing for FOF\n");
    }

    MPIU_stats(pm_comm(pm), p->np, "<->s", &npmin, &npmean, &npmax, &npstd);

    fastpm_info("load balance after first decompose : min = %g max = %g mean = %g std = %g\n",
        npmin, npmax, npmean, npstd
        );
    /* create ghosts mesh size is usually > ll so we are OK here. */
    double below[3], above[3];

    int d;
    for(d = 0; d < 3; d ++) {
        /* bigger padding reduces number of iterations */
        below[d] = -finder->linkinglength * 1;
        above[d] = finder->linkinglength * 1;
    }

    PMGhostData * pgd = pm_ghosts_create_full(pm, p,
            PACK_POS | PACK_ID | PACK_FOF, NULL,
            below, above
        );

    pm_ghosts_send(pgd, PACK_POS | PACK_ID);

    size_t np_and_ghosts = p->np + pgd->p->np;

    ptrdiff_t * head = fastpm_memory_alloc(p->mem, "FOFHead",
                    sizeof(head[0]) * np_and_ghosts, FASTPM_MEMORY_STACK);

    struct FastPMFOFData * fofsave = fastpm_memory_alloc(p->mem, "FOFSave",
                    sizeof(fofsave[0]) * np_and_ghosts, FASTPM_MEMORY_STACK);

    _fof_local_find(finder, p, pgd, head, finder->linkinglength);

    _fof_global_merge (finder, pgd, fofsave, head);

    ptrdiff_t * offset = fastpm_memory_alloc(finder->p->mem, "FOFOffset", sizeof(offset[0]) * np_and_ghosts, FASTPM_MEMORY_STACK);

    size_t nhalos = _assign_halo_attr(head, offset, np_and_ghosts, 1);

    fastpm_info("Found %td halos segments >= %d particles\n", nhalos, 1);

    fastpm_memory_free(finder->p->mem, offset);

    pm_ghosts_free(pgd);

    /* create local halos and modify head to index the local halos */
    fastpm_fof_create_local_halos(finder, halos, nhalos);

    /* now halos[head[i]] is the halo attribute of particle i. */

    /* now add the local attributes */
    fastpm_fof_compute_local_halo_attrs(finder, halos, fofsave, head, np_and_ghosts);

    fastpm_memory_free(finder->p->mem, fofsave);

    /* decompose halos by minid (gather) */
    if(0 != fastpm_store_decompose(halos,
            (fastpm_store_target_func) FastPMTargetMinID, finder, pm_comm(pm))) {

        fastpm_raise(-1, "out of space for halos decomposition.\n");
    }

    /* now head[i] is no longer the halo attribute of particle i. */

    /* local sort by minid */
    fastpm_store_sort(halos, FastPMLocalSortByMinID);

    /* reduce and update properties */
    /* set halos->mask. */
    fastpm_fof_reduce_sorted_local_halo_attrs(finder, halos);

#if 0
    if (finder->priv->ThisTask == 0)
    {
        FILE * fp;
        fp = fopen("dumpfof.data", "w");
        fwrite(halos->fof, halos->np, sizeof(halos->fof[0]), fp);
        fclose(fp);
        fp = fopen("dumpfof.mask", "w");
        fwrite(halos->mask, halos->np, sizeof(halos->mask[0]), fp);
        fclose(fp);
    }
#endif

    /* decompose halos by task (return) */
    if (0 != fastpm_store_decompose(halos,
            (fastpm_store_target_func) FastPMTargetTask, finder, pm_comm(pm))) {
        fastpm_raise(-1, "out of space for halos decomposition.\n");
    }

    /* local sort by id (restore the order) */
    fastpm_store_sort(halos, FastPMLocalSortByID);

    /* now head[i] is again the halo attribute of particle i. */

    /* the event is called with full halos, only those where mask==1 are primary
     * the others are ghosts with the correct properties but shall not show up in the
     * catalog.
     * */
    FastPMHaloEvent event[1];
    event->halos = halos;
    event->p = finder->p;
    event->ihalo = head;

    fastpm_emit_event(finder->event_handlers, FASTPM_EVENT_HALO,
                    FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event, finder);

    fastpm_memory_free(finder->p->mem, head);

    fastpm_store_subsample(halos, halos->mask, halos);
}

static void
fastpm_fof_create_local_halos(FastPMFOFFinder * finder, FastPMStore * halos, size_t nhalos)
{

    enum FastPMPackFields attributes = finder->p->attributes;
    attributes |= PACK_LENGTH | PACK_FOF;
    attributes &= ~PACK_ACC;
    attributes &= ~PACK_POTENTIAL;
    attributes &= ~PACK_DENSITY;
    attributes &= ~PACK_TIDAL;

    /* store initial position only for periodic case. non-periodic suggests light cone and
     * we cannot infer q from ID sensibly. (crashes there) */
    if(finder->priv->boxsize) {
        attributes |= PACK_Q;
    } else {
        attributes &= ~PACK_Q;
    }

    double avg_halos;
    /* + 1 to ensure avg_halos > 0 */
    MPIU_stats(pm_comm(finder->pm), nhalos + 1, "-", &avg_halos);

    /* give it enough space for rebalancing. */
    fastpm_store_init(halos, nhalos < 2 * avg_halos?2 * avg_halos:nhalos,
            attributes,
            FASTPM_MEMORY_HEAP);

    halos->np = nhalos;
    halos->a_x = finder->p->a_x;
    halos->a_v = finder->p->a_v;

    MPI_Allreduce(MPI_IN_PLACE, &nhalos, 1, MPI_LONG, MPI_SUM, pm_comm(finder->pm));
}

/* head is the hid of each particle; fofsave has the minid of each particle (unique label of each halo) */
static void
fastpm_fof_compute_local_halo_attrs (FastPMFOFFinder * finder, FastPMStore * halos, struct FastPMFOFData * fofsave, ptrdiff_t * head, size_t np_and_ghosts)
{
    FastPMStore h1[1];
    fastpm_store_init(h1, 1, halos->attributes, FASTPM_MEMORY_HEAP);

    /* */
    ptrdiff_t i;
    for(i = 0; i < halos->np; i++) {
        int d;
        halos->length[i] = 0;
        halos->id[i] = i;

        halos->mask[i] = 0; /* unselect the halos ; will turn this only if any particle is used */

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

        /* if minid is not set, this halo has no particle on this rank, and can be removed. */
        halos->fof[i].task = -1;
    }

    /* set minid and task of the halo; all of the halo particles of the same minid needs to be reduced */
    for(i = 0; i < finder->p->np; i++) {
        /* set the minid of the halo; need to take care of the ghosts too.
         * even though we do not add them to the attributes */
        ptrdiff_t hid = head[i];
        if(halos->mask[hid] == 0) {
            /* halo will be reduced by minid */
            halos->fof[hid].minid = fofsave[i].minid;
            /* halo will be returned to this task */
            halos->fof[hid].task = finder->priv->ThisTask;
            halos->mask[hid] = 1;
        } else {
            if(halos->fof[hid].minid != fofsave[i].minid) {
                abort();
            }
        }
    }

    for(i = 0; i < finder->p->np; i++) {
        ptrdiff_t hid = head[i];
        if(hid < 0) continue;

        if(halos->fof[hid].minid != fofsave[i].minid) {
            fastpm_raise(-1, "minid of halo disagrees with the particle. This shall not happend\n");
        }

        if(hid >= halos->np) {
            fastpm_raise(-1, "halo of a particle out of bounds (%td > %td)\n", hid, halos->np);
        }

        fastpm_fof_convert_particle_to_halo(finder->p, i, h1);

        fastpm_fof_add_halo_attrs(finder, halos, hid, h1, 0);
    }

    fastpm_store_destroy(h1);
}

static void
fastpm_fof_add_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos, int hid, FastPMStore * h1, int i)
{
    double * boxsize = finder->priv->boxsize;

    if(halos->aemit)
        halos->aemit[hid] += h1->aemit[i];

    int d;

    for(d = 0; d < 3; d++) {
        if(halos->v)
            halos->v[hid][d] += h1->v[i][d];
        if(halos->dx1)
            halos->dx1[hid][d] += h1->dx1[i][d];
        if(halos->dx2)
            halos->dx2[hid][d] += h1->dx2[i][d];
    }

    if(halos->x) {
        for(d = 0; d < 3; d++) {
            if(boxsize) {
                halos->x[hid][d] = periodic_add(
                    halos->x[hid][d] / halos->length[hid], halos->length[hid],
                    h1->x[i][d] / h1->length[i], h1->length[i], boxsize[d]);
            } else {
                halos->x[hid][d] += h1->x[i][d];
            }
        }
    }

    if(halos->q) {
        for(d = 0; d < 3; d ++) {
            if (boxsize) {
                halos->q[hid][d] = periodic_add(
                    halos->q[hid][d] / halos->length[hid], halos->length[hid],
                    h1->q[i][d] / h1->length[i], h1->length[i], boxsize[d]);
            } else {
                halos->q[hid][d] += h1->q[i][d];
            }
        }
    }
    /* do this after the loop because x depends on the old length. */
    halos->length[hid] += h1->length[i];
}

/* convert a particle to the first halo in the halo store */
static void
fastpm_fof_convert_particle_to_halo(FastPMStore * p, ptrdiff_t i, FastPMStore * halos)
{
    int hid = 0;

    int d;

    double q[3];

    if(halos->q)
        fastpm_store_get_q_from_id(p, p->id[i], q);

    halos->length[hid] = 1;

    for(d = 0; d < 3; d++) {
        if(halos->x)
            halos->x[hid][d] = p->x[i][d];
        if(halos->v)
            halos->v[hid][d] += p->v[i][d];
        if(halos->dx1)
            halos->dx1[hid][d] += p->dx1[i][d];
        if(halos->dx2)
            halos->dx2[hid][d] += p->dx2[i][d];
        if(halos->q) {
            halos->q[hid][d] = q[d];
        }
    }

    if(halos->aemit)
        halos->aemit[hid] = p->aemit[i];
}

static void
fastpm_fof_reduce_sorted_local_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos)
{

    ptrdiff_t i;

    ptrdiff_t first = -1;
    uint64_t lastminid = 0;
    int * ind = fastpm_memory_alloc(finder->p->mem, "HaloPermutation", sizeof(ind[0]) * halos->np, FASTPM_MEMORY_STACK);

    i = 0;
    while(i < halos->np) {
        if(halos->fof[i].task != -1) break;
        ind[i] = i;
        i++;
    }
    for(; i < halos->np + 1; i++) {
        /* we use i == halos->np to terminate the last segment */
        if(first == -1 || i == halos->np || lastminid != halos->fof[i].minid) {
            if (first >= 0) {
                /* a segment ended */
                ptrdiff_t j;
                for(j = first; j < i; j ++) {
                    ind[j] = first;
                }
            }
            if (i < halos->np) {
                /* a segment started */
                halos->mask[i] = 1;
                lastminid = halos->fof[i].minid;
            }
            /* starting a segment of the same halos */
            first = i;
        } else {
            /* inside segment, simply add the ith halo to the first of the segment */
            halos->mask[i] = 0;
            fastpm_fof_add_halo_attrs(finder, halos, first, halos, i);
        }
    }

    void * fofdata = halos->fof;
    void * id = halos->id;
    void * mask = halos->mask;

    /* use permute to replicate the first halo attr to the rest:
     *
     * we do not want to replicate
     * - fof, as fof.task is the original mpi rank of the halo
     * - id, as it is the original location of the halo on the original mpi rank. 
     * - mask, whether it is primary or not
     *   */

    halos->fof = NULL;
    halos->id = NULL;
    halos->mask = NULL;

    fastpm_store_permute(halos, ind);

    halos->fof = fofdata;
    halos->id = id;
    halos->mask = mask;

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

        /* remove halos that are shorter than nmin */
        if(halos->mask[i] && halos->length[i] < finder->nmin) {
            halos->mask[i] = 0;
        }
    }

    fastpm_memory_free(finder->p->mem, ind);
}

void
fastpm_fof_destroy(FastPMFOFFinder * finder)
{
    fastpm_destroy_event_handlers(&finder->event_handlers);
    free(finder->priv);
}


