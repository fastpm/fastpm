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
fastpm_fof_compute_halo_attrs (FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t * head);

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
FastPMTargetFOF(FastPMStore * store, ptrdiff_t i, void ** userdata)
{
    struct FastPMFOFData * fof = userdata[0];
    return fof[i].task;
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
            nmerged += _merge(&fofcomm[i], &fofsave[head[i]]);
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
        fofsave[i].minid = fofcomm[i].minid;
    }

    fastpm_memory_free(p->mem, fofcomm);
}

static void
fastpm_fof_decompose(FastPMFOFFinder * finder, FastPMStore * p, PM * pm)
{
    ptrdiff_t i;

    /* initial decompose -- reduce number of ghosts */

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

    struct FastPMFOFData * fofsave = fastpm_memory_alloc(p->mem, "FOFSave",
                    sizeof(fofsave[0]) * (p->np + pgd->p->np), FASTPM_MEMORY_STACK);

    ptrdiff_t * head = fastpm_memory_alloc(p->mem, "FOFHead",
                    sizeof(head[0]) * (p->np + pgd->p->np), FASTPM_MEMORY_STACK);

    _fof_local_find(finder, p, pgd, head, finder->linkinglength);

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

    _fof_global_merge (finder, pgd, fofsave, head);

    fastpm_memory_free(p->mem, head);

    pm_ghosts_free(pgd);


#if 0
    BigFile bf = {0};
    big_file_mpi_create(&bf, fastpm_strdup_printf("dump-fof-%d", pm->NTask), pm_comm(pm));
    write_snapshot_data(p, 1, 1, FastPMSnapshotSortByID, 0, &bf, "1", pm_comm(pm));
    big_file_mpi_close(&bf, pm_comm(pm));
#endif

    /* real decompose, move particles of the same halo to the same rank */
    {
#if 0
        /* for the assertion below only. need fofcomm to be as long as np_upper. */
        p->fof = fofcomm;
        p->attributes |= PACK_FOF;

#endif
        void * userdata[1] = {fofsave};

        if(0 != fastpm_store_decompose(p,
                    (fastpm_store_target_func) FastPMTargetFOF,
                    userdata, pm_comm(pm))) {
            fastpm_raise(-1, "out of storage space decomposing for FOF\n");
        }
#if 0
        ptrdiff_t i;
        for(i = 0; i < p->np; i ++) {
            if(p->fof[i].task != finder->priv->ThisTask) abort();
        }

        p->attributes &= ~PACK_FOF;
        p->fof = NULL;
#endif
    }

    MPIU_stats(pm_comm(pm), p->np, "<->s", &npmin, &npmean, &npmax, &npstd);

    fastpm_info("load balance after second decompose : min = %g max = %g mean = %g std = %g\n",
        npmin, npmax, npmean, npstd
        );

    fastpm_memory_free(p->mem, fofsave);
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

void
fastpm_fof_execute(FastPMFOFFinder * finder, FastPMStore * halos)
{
    fastpm_fof_decompose(finder, finder->p, finder->pm);

    KDTree tree;

    ptrdiff_t * head = fastpm_memory_alloc(finder->p->mem, "FOFHead", sizeof(head[0]) * finder->p->np, FASTPM_MEMORY_STACK);

    /* redo fof on the new decomposition -- no halo cross two ranks */
    KDNode * root = _create_kdtree(&tree, finder->kdtree_thresh, &(finder->p), 1, finder->priv->boxsize);

    kd_fof(root, finder->linkinglength, head);

    _free_kdtree(&tree, root);

    /* now do the attributes */

    fastpm_fof_compute_halo_attrs(finder, halos, head);

    FastPMHaloEvent event[1];
    event->halos = halos;
    event->p = finder->p;
    event->ihalo = head;

    fastpm_emit_event(finder->event_handlers, FASTPM_EVENT_HALO,
                    FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event, finder);

    fastpm_memory_free(finder->p->mem, head);
}

static void
fastpm_fof_compute_halo_attrs (FastPMFOFFinder * finder, FastPMStore * halos, ptrdiff_t * head)
{
    ptrdiff_t * offset = fastpm_memory_alloc(finder->p->mem, "FOFOffset", sizeof(offset[0]) * finder->p->np, FASTPM_MEMORY_STACK);

    /* */
    {
        size_t nhalos = _assign_halo_attr(head, offset, finder->p->np, finder->nmin);


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

        fastpm_info("Found %td halos >= %d particles\n", nhalos, finder->nmin);
    }

    ptrdiff_t i;
    for(i = 0; i < halos->np; i++) {
        int d;
        halos->length[i] = 0;

        if(halos->mask)
            halos->mask[i] = 1; /* select the halos for output. */

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
        halos->fof[i].minid = (uint64_t) -1;
        halos->fof[i].task = finder->priv->ThisTask;
    }

    double * boxsize = finder->priv->boxsize;

    for(i = 0; i < finder->p->np; i++) {
        ptrdiff_t hid = offset[head[i]];
        if(hid < 0) continue;

        if(hid >= halos->np) {
            abort();
        }

        if(halos->aemit)
            halos->aemit[hid] += finder->p->aemit[i];

        if(halos->fof[hid].minid > finder->p->id[i]) {
            halos->fof[hid].minid = finder->p->id[i];
        }

        int d;

        for(d = 0; d < 3; d++) {
            if(halos->v)
                halos->v[hid][d] += finder->p->v[i][d];
            if(halos->dx1)
                halos->dx1[hid][d] += finder->p->dx1[i][d];
            if(halos->dx2)
                halos->dx2[hid][d] += finder->p->dx2[i][d];
        }
        if(halos->x) {
            for(d = 0; d < 3; d++) {
                if(boxsize) {
                    halos->x[hid][d] = periodic_add(
                        halos->x[hid][d] / halos->length[hid], halos->length[hid],
                        finder->p->x[i][d], 1, boxsize[d]);
                } else {
                    halos->x[hid][d] += finder->p->x[i][d];
                }
            }
        }

        if(halos->q) {
            double q[3];
            fastpm_store_get_q_from_id(finder->p, finder->p->id[i], q);
            for(d = 0; d < 3; d ++) {
                if (boxsize) {
                    halos->q[hid][d] = periodic_add(
                        halos->q[hid][d] / halos->length[hid], halos->length[hid],
                        q[d], 1, boxsize[d]);
                } else {
                    halos->q[hid][d] += q[d];
                }
            }
        }
        /* do this after the loop because x depends on the old length. */
        halos->length[hid] += 1;
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

        if(halos->id)
            halos->id[i] = halos->fof[i].minid;
    }

    /* update head to hid, for the event. */
    for(i = 0; i < finder->p->np; i++) {
        ptrdiff_t hid = offset[head[i]];
        head[i] = hid;
    }

    fastpm_memory_free(finder->p->mem, offset);
}

void
fastpm_fof_destroy(FastPMFOFFinder * finder)
{
    fastpm_destroy_event_handlers(&finder->event_handlers);
    free(finder->priv);
}


