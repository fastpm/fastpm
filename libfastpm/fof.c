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
    MPI_Comm comm;
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

    finder->priv->comm = pm_comm(pm);
    MPI_Comm comm = finder->priv->comm;
    MPI_Comm_rank(comm, &finder->priv->ThisTask);
    MPI_Comm_size(comm, &finder->priv->NTask);
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
_merge(FastPMStore * remote, ptrdiff_t i, FastPMStore * fofhead, ptrdiff_t j)
{
    int merge = 0;
    if(remote->minid[i] < fofhead->minid[j]) {
        merge = 1;
    }
    if(remote->minid[i] == fofhead->minid[j] &&
       remote->task[i] != fofhead->task[j] &&
       remote->task[i] >= 0) {
        merge = 1;
    }

    if(merge) {
        fofhead->minid[j] = remote->minid[i];
        fofhead->task[j] = remote->task[i];

    }
    return merge;
}

static void
_fof_global_merge(
    FastPMFOFFinder * finder,
    PMGhostData * pgd,
    FastPMStore * savebuff,
    ptrdiff_t * head
)
{
    ptrdiff_t i;
    FastPMStore * p = finder->p;
    MPI_Comm comm = finder->priv->comm;

    FastPMStore commbuff[1];

    fastpm_store_init(commbuff, p->np + pgd->p->np, PACK_MINID | PACK_TASK, FASTPM_MEMORY_STACK);

    size_t npmax = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &npmax, 1, MPI_LONG, MPI_MAX, comm);

    /* initialize minid, used as a global tag of groups as we merge */
    for(i = 0; i < p->np; i ++) {
        /* assign unique ID to each particle; could use a better scheme with true offsets */
        savebuff->minid[i] = i + finder->priv->ThisTask * npmax;
        savebuff->task[i] = finder->priv->ThisTask;
    }

    /* send minid */
    p->minid = savebuff->minid;
    p->task = savebuff->task;
    pm_ghosts_send(pgd, PACK_MINID | PACK_TASK);
    p->minid = NULL;
    p->task = NULL;

    for(i = 0; i < pgd->p->np; i ++) {
        savebuff->minid[i + p->np] = pgd->p->minid[i];
        savebuff->task[i + p->np] = -1; /* ghosts we do not yet know their host task */
    }

    /* reduce the minid of the head items according to the local connection. */

    for(i = 0; i < p->np + pgd->p->np; i ++) {
        _merge(savebuff, i, savebuff, head[i]);
    }

    /* at this point all items with head[i] = i have local minid and task */

    int iter = 0;

    while(1) {

        /* prepare the communication buffer, every particle has
         * the minid and task of the current head. such that
         * they will connect to the other ranks correctly */
        for(i = 0; i < p->np + pgd->p->np; i ++) {
            commbuff->minid[i] = savebuff->minid[head[i]];
            commbuff->task[i] = savebuff->task[head[i]];
        }

        /* at this point all items on fofcomm have local minid and task, ready to send */

        p->minid = commbuff->minid;
        p->task = commbuff->task;
        pm_ghosts_send(pgd, PACK_MINID | PACK_TASK);

        p->minid = NULL;
        p->task = NULL;

        /* flatten the fof data received from the ghosts */
        /* FIXME: use reduce? */ 
        for(i = p->np; i < p->np + pgd->p->np; i ++) {
            commbuff->minid[i] = pgd->p->minid[i - p->np];
            commbuff->task[i] = pgd->p->task[i - p->np];
        }
        size_t nmerged = 0;
        for(i = p->np; i < p->np + pgd->p->np; i ++) {
            int m = _merge(commbuff, i, savebuff, head[i]);
            nmerged += m;
        }

        /* at this point all items on savebuff->fof with head[i] = i have present minid and task */

        MPI_Allreduce(MPI_IN_PLACE, &nmerged, 1, MPI_LONG, MPI_SUM, comm);

        MPI_Barrier(comm);

        fastpm_info("FOF reduction iteration %d : merged %td crosslinks\n", iter, nmerged);

        if(nmerged == 0) break;

        for(i = 0; i < p->np + pgd->p->np; i ++) {
            if(savebuff->minid[i] < savebuff->minid[head[i]]) {
                fastpm_raise(-1, "savebuff->fof invariance is broken i = %td np = %td\n", i, p->np);
            }
        }

        iter++;
    }

    for(i = 0; i < p->np; i ++) {
        if(commbuff->task[i] == -1) {
            fastpm_raise(-1, "undeterimined particle %d id = %03td head = %03td head[%03d/%03d] : id = %03td task = %03d\n",
                finder->priv->ThisTask, p->id[i], p->id[head[i]], i, p->np, commbuff->minid[i], commbuff->task[i]);
        }
        savebuff->minid[i] = commbuff->minid[i];
        savebuff->task[i] = commbuff->task[i];
    }
    for(i = p->np; i < p->np + pgd->p->np; i ++) {
        savebuff->task[i] = -1;
    }

    fastpm_store_destroy(commbuff);
}

/* set head[i] to hid*/
static size_t
_assign_halo_attr(FastPMFOFFinder * finder, ptrdiff_t * head, size_t np, size_t np_ghosts, int nmin)
{
    ptrdiff_t * offset = fastpm_memory_alloc(finder->p->mem, "FOFOffset", sizeof(offset[0]) * (np + np_ghosts), FASTPM_MEMORY_STACK);
    uint8_t * has_remote = fastpm_memory_alloc(finder->p->mem, "FOFHasRemote", sizeof(has_remote[0]) * (np + np_ghosts), FASTPM_MEMORY_STACK);
    uint8_t * has_local = fastpm_memory_alloc(finder->p->mem, "FOFHasLocal", sizeof(has_local[0]) * (np + np_ghosts), FASTPM_MEMORY_STACK);

    ptrdiff_t i;
    for(i = 0; i < np + np_ghosts; i ++) {
        offset[i] = 0;
        has_remote[i] = 0;
        has_local[i] = 0;
    }

    /* set offset to number of particles in the halo */
    for(i = 0; i < np + np_ghosts; i ++) {
        offset[head[i]] ++;
    }

    for(i = 0; i < np; i ++) {
        has_local[head[i]] = 1;
    }
    /* if the group is connected to a remote component */
    for(i = np; i < np + np_ghosts; i ++) {
        has_remote[head[i]] = 1;
    }

    size_t it = 0;

    /* assign attr index for groups at least contain 1 local particle */
    for(i = 0; i < np + np_ghosts; i ++) {
        if(has_local[i] && (has_remote[i] || offset[i] >= nmin)) {
            offset[i] = it;
            it ++;
        } else {
            offset[i] = -1;
        }
    }

    size_t nhalos = it;

    /* update head [i] to offset[head[i]], which stores the index in the halo store for this particle. */
    for(i = 0; i < np + np_ghosts; i ++) {
        head[i] = offset[head[i]];
        /* this will not happen if nmin == 1 */
        if(head[i] != -1 && head[i] > nhalos) {
            fastpm_raise(-1, "head[i] (%td) > nhalos (%td) This shall not happen.\n", head[i], nhalos);
        }
    }

    fastpm_memory_free(finder->p->mem, has_local);
    fastpm_memory_free(finder->p->mem, has_remote);
    fastpm_memory_free(finder->p->mem, offset);

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

/*
 * apply mask to a halo storage.
 * relable head, and halo->id
 * */
static void
fastpm_fof_subsample(FastPMFOFFinder * finder, FastPMStore * halos, uint8_t * mask, ptrdiff_t * head)
{
    /* mapping goes from the old head[i] value to the new head[i] value */
    ptrdiff_t * mapping = fastpm_memory_alloc(finder->p->mem, "FOFMapping", sizeof(mapping[0]) * halos->np, FASTPM_MEMORY_STACK);

    ptrdiff_t i;
    for(i = 0; i < halos->np; i ++) {
        mapping[i] = -1;
        halos->id[i] = i;
    }

    /* remove non-contributing halo segments */
    fastpm_store_subsample(halos, mask, halos);

    for(i = 0; i < halos->np; i ++) {
        mapping[halos->id[i]] = i;
        halos->id[i] = i;
    }

    /* adjust head[i] */

    for(i = 0; i < finder->p->np; i ++) {
        if(head[i] >= 0) {
            head[i] = mapping[head[i]];
        }
    }

    fastpm_memory_free(finder->p->mem, mapping);
}

static void
fastpm_fof_apply_length_cut(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t * head)
{
    MPI_Comm comm = finder->priv->comm;

    uint8_t * mask = fastpm_memory_alloc(finder->p->mem, "LengthMask", sizeof(mask[0]) * halos->np, FASTPM_MEMORY_STACK);
    ptrdiff_t i;
    for(i = 0; i < halos->np; i ++) {
        /* remove halos that are shorter than nmin */
        if(halos->length[i] < finder->nmin) {
            mask[i] = 0;
        } else {
            mask[i] = 1;
        }
    }
    fastpm_fof_subsample(finder, halos, mask, head);
    fastpm_memory_free(finder->p->mem, mask);

    fastpm_info("After length cut we have %td halos (including ghost halos)\n", fastpm_store_get_np_total(halos, comm));
}

/* first every undecided halo; then the real halo with particles */
static int
FastPMLocalSortByMinID(const int i1,
                    const int i2,
                    FastPMStore * p)
{
    int v1 = (p->minid[i1] < p->minid[i2]);
    int v2 = (p->minid[i1] > p->minid[i2]);

    return v2 - v1;
}

static int
FastPMTargetMinID(FastPMStore * store, ptrdiff_t i, void * userdata)
{
    FastPMFOFFinder * finder = userdata;

    const uint32_t GOLDEN32 = 2654435761ul;
    /* const uint64_t GOLDEN64 = 11400714819323198549; */
    /* may over flow, but should be okay here as the periodicity is ggt NTask */
    int key = (store->minid[i] * GOLDEN32) % (unsigned) finder->priv->NTask;
    return key;
}

static int
FastPMTargetTask(FastPMStore * store, ptrdiff_t i, void * userdata)
{
    return store->task[i];
}

/*
 * This function group halos by halos->minid[i].
 *
 * All halos with the same minid will have the same attribute values afterwards, except
 * a few book keeping items named in this function (see comments inside)
 *
 * */
static void
fastpm_fof_reduce_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos,
        void (* add_func) (FastPMFOFFinder * finder, FastPMStore * halo1, ptrdiff_t i1, FastPMStore * halo2, ptrdiff_t i2 ),
        void (* reduce_func)(FastPMFOFFinder * finder, FastPMStore * halo1, ptrdiff_t i1)
)
{

    ptrdiff_t i;

    ptrdiff_t first = -1;
    uint64_t lastminid = 0;

    /* ind is the array to use to replicate items */
    int * ind = fastpm_memory_alloc(finder->p->mem, "HaloPermutation", sizeof(ind[0]) * halos->np, FASTPM_MEMORY_STACK);

    /* the following items will have mask[i] == 1, but we will mark some to 0 if
     * they are not the principle (first) halos segment with this minid */
    for(i = 0; i < halos->np + 1; i++) {
        /* we use i == halos->np to terminate the last segment */
        if(first == -1 || i == halos->np || lastminid != halos->minid[i]) {
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
                lastminid = halos->minid[i];
            }
            /* starting a segment of the same halos */
            first = i;
        } else {
            /* inside segment, simply add the ith halo to the first of the segment */
            halos->mask[i] = 0;
            add_func(finder, halos, first, halos, i);
        }
    }

    void * minid = halos->minid;
    void * task = halos->task;
    void * id = halos->id;
    void * mask = halos->mask;

    /* use permute to replicate the first halo attr to the rest:
     *
     * we do not want to replicate
     * - fof, as fof.task is the original mpi rank of the halo
     * - id, as it is the original location of the halo on the original mpi rank. 
     * - mask, whether it is primary or not
     *
     * we need to replicate because otherwise when we return the head array on the
     * original ranks will be violated.
     *   */

    halos->minid = NULL;
    halos->task = NULL;
    halos->id = NULL;
    halos->mask = NULL;

    fastpm_store_permute(halos, ind);

    halos->minid = minid;
    halos->task = task;
    halos->id = id;
    halos->mask = mask;

    for(i = 0; i < halos->np; i++) {
        reduce_func(finder, halos, i);
    }

    fastpm_memory_free(finder->p->mem, ind);
}


static void
fastpm_fof_compute_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos,
            ptrdiff_t * head,
            void (*convert_func)(FastPMFOFFinder * finder, FastPMStore * p, ptrdiff_t i, FastPMStore * halos),
            void (*add_func)(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t i1, FastPMStore * halos2, ptrdiff_t i2),
            void (*reduce_func)(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t i)
)
{
    MPI_Comm comm = finder->priv->comm;

    FastPMStore h1[1];
    fastpm_store_init(h1, 1, halos->attributes, FASTPM_MEMORY_HEAP);
    ptrdiff_t i;

    for(i = 0; i < finder->p->np; i++) {
        ptrdiff_t hid = head[i];
        if(hid < 0) continue;


        if(hid >= halos->np) {
            fastpm_raise(-1, "halo of a particle out of bounds (%td > %td)\n", hid, halos->np);
        }

        /* initialize h1 with existing halo attributes for this particle */
        fastpm_store_take(halos, hid, h1, 0);

        convert_func(finder, finder->p, i, h1);

        add_func(finder, halos, hid, h1, 0);
    }

    fastpm_store_destroy(h1);
    /* decompose halos by minid (gather); if all halo segments of the same minid are on the same rank, we can combine
     * these into a single entry, then replicate for each particle to look up;
     * halo segments that have no local particles are never exchanged to another rank. */
    if(0 != fastpm_store_decompose(halos,
            (fastpm_store_target_func) FastPMTargetMinID, finder, comm)) {

        fastpm_raise(-1, "out of space sending halos by MinID.\n");
    }

    /* now head[i] is no longer the halo attribute of particle i. */

    /* to combine, first local sort by minid; those without local particles are moved to the beginning
     * so we can easily skip them. */
    fastpm_store_sort(halos, FastPMLocalSortByMinID);

    /* reduce and update properties */
    fastpm_fof_reduce_halo_attrs(finder, halos, add_func, reduce_func);

    /* decompose halos by task (return) */
    if (0 != fastpm_store_decompose(halos,
            (fastpm_store_target_func) FastPMTargetTask, finder, comm)) {
        fastpm_raise(-1, "out of space for gathering halos this shall never happen.\n");
    }
    /* local sort by id (restore the order) */
    fastpm_store_sort(halos, FastPMLocalSortByID);

    /* now head[i] is again the halo attribute of particle i. */
}

/*
 * compute the attrs of the local halo segments based on local particles.
 * head : the hid of each particle; 
 * fofsave : the minid of each particle (unique label of each halo)
 *
 * if a halo has no local particles,  halos->mask[i] is set to 0.
 * if a halo has any local particles, and halos->mask[i] is set to 1.
 * */
static void
fastpm_fof_remove_empty_halos(FastPMFOFFinder * finder, FastPMStore * halos, FastPMStore * savebuff, ptrdiff_t * head)
{

    /* */
    ptrdiff_t i;

    /* set minid and task of the halo; all of the halo particles of the same minid needs to be reduced */
    for(i = 0; i < finder->p->np; i++) {
        /* set the minid of the halo; need to take care of the ghosts too.
         * even though we do not add them to the attributes */
        ptrdiff_t hid = head[i];

        if(hid < 0) continue;

        if(halos->mask[hid] == 0) {
            /* halo will be reduced by minid */
            halos->minid[hid] = savebuff->minid[i];
            halos->mask[hid] = 1;
        } else {
            if(halos->minid[hid] != savebuff->minid[i]) {
                fastpm_raise(-1, "Consistency check failed after FOF global merge.\n");
            }
        }
    }

    /* halo will be returned to this task;
     * we save ThisTask here.*/
    for(i = 0; i < halos->np; i ++) {
        halos->task[i] = finder->priv->ThisTask;
    }

    fastpm_fof_subsample(finder, halos, halos->mask, head);
}

static void
_add_basic_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t hid, FastPMStore * h1, ptrdiff_t i)
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
_convert_basic_halo_attrs(FastPMFOFFinder * finder, FastPMStore * p, ptrdiff_t i, FastPMStore * halos)
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
            halos->v[hid][d] = p->v[i][d];
        if(halos->dx1)
            halos->dx1[hid][d] = p->dx1[i][d];
        if(halos->dx2)
            halos->dx2[hid][d] = p->dx2[i][d];
        if(halos->q) {
            halos->q[hid][d] = q[d];
        }
    }

    if(halos->aemit)
        halos->aemit[hid] = p->aemit[i];
}
static void
_reduce_basic_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t i)
{
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

static void
_add_extended_halo_attrs(FastPMFOFFinder * finder, FastPMStore * h1, ptrdiff_t i1, FastPMStore * h2, ptrdiff_t i2)
{
    int d;
    if(h1->rvdisp) {
        for(d = 0; d < 9; d ++) {
            h1->rvdisp[i1][d] += h2->rvdisp[i2][d];
        }
    }
    if(h1->vdisp) {
        for(d = 0; d < 6; d ++) {
            h1->vdisp[i1][d] += h2->vdisp[i2][d];
        }
    }
    if(h1->rdisp) {
        for(d = 0; d < 6; d ++) {
            h1->rdisp[i1][d] += h2->rdisp[i2][d];
        }
    }
}

static void
_convert_extended_halo_attrs(FastPMFOFFinder * finder, FastPMStore * p, ptrdiff_t i, FastPMStore * halos)
{
    int hid = 0;

    int d;
    double rrel[3];

    for(d = 0; d < 3; d ++) {
        rrel[d] = p->x[i][d] - halos->x[hid][d];

        if(finder->priv->boxsize) {
            double L = finder->priv->boxsize[d];
            while(rrel[d] > L / 2) rrel[d] -= L;
            while(rrel[d] < -L / 2) rrel[d] += L;
        }
    }

    if(halos->vdisp) {
        /* FIXME: add hubble expansion term based on aemit and hubble function? needs to modify Finder object */

        double vrel[3];
        for(d = 0; d < 3; d ++) {
            vrel[d] = p->v[i][d] - halos->v[hid][d];
        }
        for(d = 0; d < 3; d ++) {
            halos->vdisp[hid][d] = vrel[d] * vrel[d];
            halos->vdisp[hid][d + 3] = vrel[d] * vrel[(d + 1) % 3];
        }
    }
    if(halos->rvdisp) {
        /* FIXME: add hubble expansion term based on aemit and hubble function? needs to modify Finder object */

        double vrel[3];
        for(d = 0; d < 3; d ++) {
            vrel[d] = p->v[i][d] - halos->v[hid][d];
        }
        for(d = 0; d < 3; d ++) {
            halos->rvdisp[hid][d] = rrel[d] * vrel[d];
            halos->rvdisp[hid][d + 3] = rrel[d] * vrel[(d + 1) % 3];
            halos->rvdisp[hid][d + 6] = rrel[d] * vrel[(d + 2) % 3];
        }
    }
    if(halos->rdisp) {
        for(d = 0; d < 3; d ++) {
            halos->rdisp[hid][d] = rrel[d] * rrel[d];
            halos->rdisp[hid][d + 3] = rrel[d] * rrel[(d + 1) % 3];
        }
    }
}
static void
_reduce_extended_halo_attrs(FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t hid)
{
    double n = halos->length[hid];
    int d;
    if(halos->rvdisp) {
        for(d = 0; d < 9; d ++) {
            halos->rvdisp[hid][d] /= n;
        }
    }
    if(halos->vdisp) {
        for(d = 0; d < 6; d ++) {
            halos->vdisp[hid][d] /= n;
        }
    }
    if(halos->rdisp) {
        for(d = 0; d < 3; d ++) {
            halos->rdisp[hid][d] /= n;
        }
    }
}

/* This function creates the storage object for halo segments that are local on this
 * rank. We store mnay attributes. We only allow a flucutation of 2 around avg_halos.
 * this should be OK, since we will only redistribute by the MinID, which are supposed
 * to be very uniform.
 * */
static void
fastpm_fof_create_local_halos(FastPMFOFFinder * finder, FastPMStore * halos, size_t nhalos)
{

    MPI_Comm comm = finder->priv->comm;

    enum FastPMPackFields attributes = finder->p->attributes;
    attributes |= PACK_LENGTH | PACK_MINID | PACK_TASK;
    attributes |= PACK_RDISP | PACK_VDISP | PACK_RVDISP;
    attributes |= PACK_ACC; /* ACC used as the first particle position offset */
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
    double max_halos;
    /* + 1 to ensure avg_halos > 0 */
    MPIU_stats(comm, nhalos + 1, "->", &avg_halos, &max_halos);

    /* give it enough space for rebalancing. */
    fastpm_store_init(halos, (size_t) (max_halos * 2),
            attributes,
            FASTPM_MEMORY_HEAP);

    halos->np = nhalos;
    halos->a_x = finder->p->a_x;
    halos->a_v = finder->p->a_v;

    ptrdiff_t i;
    for(i = 0; i < halos->np; i++) {
        halos->id[i] = i;

        halos->mask[i] = 0; /* unselect the halos ; will turn this only if any particle is used */
        /* everthing should have been set to zero already by fastpm_store_init */
    }
}


void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos)
{
    /* initial decompose -- reduce number of ghosts */
    FastPMStore * p = finder->p;
    PM * pm = finder->pm;
    MPI_Comm comm = finder->priv->comm;

    /* only do wrapping for periodic data */
    if(finder->priv->boxsize)
        fastpm_store_wrap(p, finder->priv->boxsize);

    double npmax, npmin, npstd, npmean;

    MPIU_stats(comm, p->np, "<->s", &npmin, &npmean, &npmax, &npstd);

    fastpm_info("load balance after decompose : min = %g max = %g mean = %g std = %g\n",
        npmin, npmax, npmean, npstd
        );

    /* still route particles to the pm pencils as if they are periodic. */
    if(0 != fastpm_store_decompose(p,
                (fastpm_store_target_func) FastPMTargetPM,
                pm, comm)
    ) {
        fastpm_raise(-1, "out of storage space decomposing for FOF\n");
    }

    MPIU_stats(comm, p->np, "<->s", &npmin, &npmean, &npmax, &npstd);

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
            PACK_POS | PACK_ID | PACK_MINID | PACK_TASK, NULL,
            below, above
        );

    pm_ghosts_send(pgd, PACK_POS | PACK_ID);

    size_t np_and_ghosts = p->np + pgd->p->np;

    ptrdiff_t * head = fastpm_memory_alloc(p->mem, "FOFHead",
                    sizeof(head[0]) * np_and_ghosts, FASTPM_MEMORY_STACK);

    FastPMStore savebuff[1];
    fastpm_store_init(savebuff, np_and_ghosts, PACK_MINID | PACK_TASK, FASTPM_MEMORY_STACK);

    _fof_local_find(finder, p, pgd, head, finder->linkinglength);

    _fof_global_merge (finder, pgd, savebuff, head);

    /* assign halo attr entries. This will keep only candidates that can possibly reach to nmin */
    size_t nhalos = _assign_halo_attr(finder, head, p->np, pgd->p->np, finder->nmin);

    fastpm_info("Found %td halos segments >= %d particles; or cross linked. \n", nhalos, finder->nmin);

    pm_ghosts_free(pgd);

    /* create local halos and modify head to index the local halos */
    fastpm_fof_create_local_halos(finder, halos, nhalos);

    /* remove halos without any local particles */
    fastpm_fof_remove_empty_halos(finder, halos, savebuff, head);

    fastpm_store_destroy(savebuff);

    /* reduce the primary halo attrs */
    fastpm_fof_compute_halo_attrs(finder, halos, head, _convert_basic_halo_attrs, _add_basic_halo_attrs, _reduce_basic_halo_attrs);

    /* apply length cut */
    fastpm_fof_apply_length_cut(finder, halos, head);

    /* reduce the primary halo attrs */
    fastpm_fof_compute_halo_attrs(finder, halos, head, _convert_extended_halo_attrs, _add_extended_halo_attrs, _reduce_extended_halo_attrs);

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

    fastpm_info("After event: %td halos.\n", fastpm_store_get_np_total(halos, comm));
}

void
fastpm_fof_destroy(FastPMFOFFinder * finder)
{
    fastpm_destroy_event_handlers(&finder->event_handlers);
    free(finder->priv);
}


