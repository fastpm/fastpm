#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <alloca.h>
#include <string.h>
#include "bigfile-mpi.h"
#include "bigfile-internal.h"

/* disable aggregation by default */
static size_t _BigFileAggThreshold = 0;
static int _big_file_mpi_verbose = 0;

static int big_block_mpi_broadcast(BigBlock * bb, int root, MPI_Comm comm);
static int big_file_mpi_broadcast_anyerror(int rt, MPI_Comm comm);

#define BCAST_AND_RAISEIF(rt, comm) \
    if(0 != (rt = big_file_mpi_broadcast_anyerror(rt, comm))) { \
        return rt; \
    } \

void
big_file_mpi_set_verbose(int verbose)
{
    _big_file_mpi_verbose = verbose;
}

void
big_file_mpi_set_aggregated_threshold(size_t bytes)
{
    _BigFileAggThreshold = bytes;
}

size_t
big_file_mpi_get_aggregated_threshold()
{
    return _BigFileAggThreshold;
}

int big_file_mpi_open(BigFile * bf, const char * basename, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int rt = 0;
    if (rank == 0) {
        rt = big_file_open(bf, basename);
    } else {
        /* FIXME : */
        bf->basename = _strdup(basename);
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);

    return rt;
}

int big_file_mpi_create(BigFile * bf, const char * basename, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int rt;
    if (rank == 0) {
        rt = big_file_create(bf, basename);
    } else {
        /* FIXME : */
        bf->basename = _strdup(basename);
        rt = 0;
    }
    BCAST_AND_RAISEIF(rt, comm);

    return rt;
}

/**Helper function for big_file_mpi_create_block, above*/
static int _big_block_mpi_create(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, size_t fsize[], MPI_Comm comm);

/** Helper function for big_file_mpi_open_block, above*/
static int _big_block_mpi_open(BigBlock * bb, const char * basename, MPI_Comm comm);

int big_file_mpi_open_block(BigFile * bf, BigBlock * block, const char * blockname, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    if(!bf || !bf->basename || !blockname) return 1;
    char * basename = alloca(strlen(bf->basename) + strlen(blockname) + 128);
    sprintf(basename, "%s/%s/", bf->basename, blockname);
    return _big_block_mpi_open(block, basename, comm);
}

int big_file_mpi_create_block(BigFile * bf, BigBlock * block, const char * blockname, const char * dtype, int nmemb, int Nfile, size_t size, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;

    size_t fsize[Nfile];
    int i;
    for(i = 0; i < Nfile; i ++) {
        fsize[i] = size * (i + 1) / Nfile 
                 - size * (i) / Nfile;
    }
    int rank;
    MPI_Comm_rank(comm, &rank);

    int rt = 0;
    if (rank == 0) {
        rt = _big_file_mksubdir_r(bf->basename, blockname);
    } else {
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);

    char * basename = alloca(strlen(bf->basename) + strlen(blockname) + 128);
    sprintf(basename, "%s/%s/", bf->basename, blockname);
    return _big_block_mpi_create(block, basename, dtype, nmemb, Nfile, fsize, comm);
}

int big_file_mpi_close(BigFile * bf, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rt = big_file_close(bf);
    MPI_Barrier(comm);
    return rt;
}

static int
_big_block_mpi_open(BigBlock * bb, const char * basename, MPI_Comm comm)
{
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int rt;
    if(rank == 0) { 
        rt = _big_block_open(bb, basename);
    } else {
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);

    big_block_mpi_broadcast(bb, 0, comm);
    return 0;
}

static int
_big_block_mpi_create(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, size_t fsize[], MPI_Comm comm)
{
    int rank;
    int NTask;
    int rt;

    if(comm == MPI_COMM_NULL) return 0;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &rank);

    if(rank == 0) {
        rt = _big_block_create_internal(bb, basename, dtype, nmemb, Nfile, fsize);
    } else {
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);

    big_block_mpi_broadcast(bb, 0, comm);

    int i;
    for(i = (size_t) bb->Nfile * rank / NTask; i < (size_t) bb->Nfile * (rank + 1) / NTask; i ++) {
        FILE * fp = _big_file_open_a_file(bb->basename, i, "w", 1);
        if(fp == NULL) {
            rt = -1;
            break;
        }
        fclose(fp);
    }

    BCAST_AND_RAISEIF(rt, comm);

    return rt;
}

int big_block_mpi_grow(BigBlock * bb, int Nfile_grow, size_t fsize_grow[], MPI_Comm comm) {
    int rank;
    int NTask;
    int rt;

    if(comm == MPI_COMM_NULL) return 0;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &rank);

    int oldNfile = bb->Nfile;

    if(rank == 0) {
        rt = _big_block_grow_internal(bb, Nfile_grow, fsize_grow);
    } else {
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);

    if(rank != 0) {
        /* closed on non-root because we will bcast.*/
        _big_block_close_internal(bb);
    }
    big_block_mpi_broadcast(bb, 0, comm);

    int i;
    for(i = (size_t) Nfile_grow * rank / NTask; i < (size_t) Nfile_grow * (rank + 1) / NTask; i ++) {
        FILE * fp = _big_file_open_a_file(bb->basename, i + oldNfile, "w", 1);
        if(fp == NULL) {
            rt = -1;
            break;
        }
        fclose(fp);
    }

    BCAST_AND_RAISEIF(rt, comm);

    return rt;
}

int
big_block_mpi_grow_simple(BigBlock * bb, int Nfile_grow, size_t size_grow, MPI_Comm comm)
{
    size_t fsize[Nfile_grow];
    int i;
    for(i = 0; i < Nfile_grow; i ++) {
        fsize[i] = size_grow * (i + 1) / Nfile_grow
                 - size_grow * (i) / Nfile_grow;
    }
    int rank;
    MPI_Comm_rank(comm, &rank);

    return big_block_mpi_grow(bb, Nfile_grow, fsize, comm);
}


int
big_block_mpi_flush(BigBlock * block, MPI_Comm comm)
{
    if(comm == MPI_COMM_NULL) return 0;

    int rank;
    MPI_Comm_rank(comm, &rank);

    unsigned int * checksum = alloca(sizeof(int) * block->Nfile);
    MPI_Reduce(block->fchecksum, checksum, block->Nfile, MPI_UNSIGNED, MPI_SUM, 0, comm);
    int dirty;
    MPI_Reduce(&block->dirty, &dirty, 1, MPI_INT, MPI_LOR, 0, comm);
    int rt;
    if(rank == 0) {
        /* only the root rank updates */
        int i;
        big_block_set_dirty(block, dirty);
        for(i = 0; i < block->Nfile; i ++) {
            block->fchecksum[i] = checksum[i];
        }
        rt = big_block_flush(block);
    } else {
        rt = 0;
    }

    BCAST_AND_RAISEIF(rt, comm);
    /* close as we will broadcast the block */
    if(rank != 0) {
        _big_block_close_internal(block);
    }
    big_block_mpi_broadcast(block, 0, comm);
    return 0;

}
int big_block_mpi_close(BigBlock * block, MPI_Comm comm) {

    int rt = big_block_mpi_flush(block, comm);
    _big_block_close_internal(block);

    return rt;
}

static int
big_file_mpi_broadcast_anyerror(int rt, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    struct {
        int value;
        int loc;
    } ii = {rt != 0, rank};

    MPI_Allreduce(MPI_IN_PLACE, &ii, 1, MPI_2INT, MPI_MAXLOC, comm);

    if (ii.value == 0) {
        /* no errors */
        return 0;
    }
    int root = ii.loc;

    char * error = big_file_get_error_message();

    int errorlen;
    if(rank == root) {
        errorlen = strlen(error);
    }
    MPI_Bcast(&errorlen, 1, MPI_INT, root, comm);

    if(rank != root) {
        error = malloc(errorlen + 1);
    }

    MPI_Bcast(error, errorlen + 1, MPI_BYTE, root, comm);

    if(rank != root) {
        big_file_set_error_message(error);
        free(error);
    }

    MPI_Bcast(&rt, 1, MPI_INT, root, comm);

    return rt;
}

static int
big_block_mpi_broadcast(BigBlock * bb, int root, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    void * buf = NULL;
    size_t bytes = 0;

    if(rank == root) {
        buf = _big_block_pack(bb, &bytes);
    }

    MPI_Bcast(&bytes, sizeof(bytes), MPI_BYTE, root, comm);

    if(rank != root) {
        buf = malloc(bytes);
    }

    MPI_Bcast(buf, bytes, MPI_BYTE, root, comm);

    _big_block_unpack(bb, buf);

    free(buf);
    return 0;
}

static int
_assign_colors(size_t glocalsize, size_t * sizes, int * ncolor, MPI_Comm comm)
{
    int NTask;
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    int i;
    int mycolor;
    size_t current_size = 0;
    int current_color = 0;
    int lastcolor = 0;
    for(i = 0; i < NTask; i ++) {
        lastcolor = current_color;

        if(i == ThisTask) {
            mycolor = lastcolor;
        }

        if (_big_file_mpi_verbose) {
            if(ThisTask == 0) {
                printf("current_size %td, sizes[i] %td, glocalsize %td, current_color %d\n",
                    current_size, sizes[i], glocalsize, current_color);
            }
        }
        current_size += sizes[i];

        if(current_size >= glocalsize) {
            current_size = 0;
            current_color ++;
        }
    }
    /* no data for color of -1; exclude them later with special cases */
    if(sizes[ThisTask] == 0) {
        mycolor = -1;
    }
    *ncolor = lastcolor + 1;
    return mycolor;
}

static size_t
_collect_sizes(size_t localsize, size_t * sizes, size_t * myoffset, MPI_Comm comm)
{

    int ThisTask, NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    size_t totalsize;

    sizes[ThisTask] = localsize;

    MPI_Datatype MPI_PTRDIFFT;

    if(sizeof(ptrdiff_t) == sizeof(long)) {
        MPI_PTRDIFFT = MPI_LONG;
    } else if(sizeof(ptrdiff_t) == sizeof(int)) {
        MPI_PTRDIFFT = MPI_INT;
    } else { abort(); }

    MPI_Allreduce(&sizes[ThisTask], &totalsize, 1, MPI_PTRDIFFT, MPI_SUM, comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sizes, 1, MPI_PTRDIFFT, comm);

    int i;
    *myoffset = 0;
    for(i = 0; i < ThisTask; i ++) {
        (*myoffset) += sizes[i];
    }

    return totalsize;
}

struct SegmentGroupDescr {
    /* data model: rank <- segment <- group */
    int Ngroup;
    int Nsegments;
    int GroupID; /* ID of the group of this rank */
    int ThisSegment; /* SegmentID of the local data chunk on this rank*/

    int segment_start; /* segments responsible in this group */
    int segment_end;

    int is_group_leader;
    int group_leader_rank;
    int segment_leader_rank;
    MPI_Comm Group;  /* communicator for all ranks in the group */
    MPI_Comm Leader; /* communicator for all ranks that are group leaders */
    MPI_Comm Segment; /* communicator for all ranks in this segment */
};

static void
_create_segment_group(struct SegmentGroupDescr * descr, size_t * sizes, size_t avgsegsize, int Ngroup, MPI_Comm comm)
{
    int i;
    int ThisTask, NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    descr->ThisSegment = _assign_colors(avgsegsize, sizes, &descr->Nsegments, comm);

    if(descr->ThisSegment >= 0) {
        /* assign segments to groups.
         * if Nsegments < Ngroup, some groups will have no segments, and thus no ranks belong to them. */
        descr->GroupID = ((size_t) descr->ThisSegment) * Ngroup / descr->Nsegments;
    } else {
        descr->GroupID = Ngroup + 1;
        descr->ThisSegment = NTask + 1;
    }

    descr->Ngroup = Ngroup;

    MPI_Comm_split(comm, descr->GroupID, ThisTask, &descr->Group);

    MPI_Allreduce(&descr->ThisSegment, &descr->segment_start, 1, MPI_INT, MPI_MIN, descr->Group);
    MPI_Allreduce(&descr->ThisSegment, &descr->segment_end, 1, MPI_INT, MPI_MAX, descr->Group);

    descr->segment_end ++;

    int rank;

    MPI_Comm_rank(descr->Group, &rank);

    struct {
        size_t val;
        int   rank;
    } leader_st;

    leader_st.val = sizes[ThisTask];
    leader_st.rank = rank;

    MPI_Allreduce(MPI_IN_PLACE, &leader_st, 1, MPI_LONG_INT, MPI_MAXLOC, descr->Group);

    descr->is_group_leader = rank == leader_st.rank;
    descr->group_leader_rank = leader_st.rank;

    MPI_Comm_split(comm, rank == leader_st.rank? 0 : 1, ThisTask, &descr->Leader);

    MPI_Comm_split(descr->Group, descr->ThisSegment, ThisTask, &descr->Segment);
    int rank2;

    MPI_Comm_rank(descr->Segment, &rank2);

    leader_st.val = sizes[ThisTask];
    leader_st.rank = rank2;

    MPI_Allreduce(MPI_IN_PLACE, &leader_st, 1, MPI_LONG_INT, MPI_MINLOC, descr->Segment);
    descr->segment_leader_rank = leader_st.rank;

    if (_big_file_mpi_verbose) {
        printf("bigfile: ThisTask = %d ThisSegment = %d / %d GroupID = %d / %d rank = %d rank in segment = %d segstart = %d segend = %d segroot = %d\n",
            ThisTask, descr->ThisSegment, descr->Nsegments, descr->GroupID, descr->Ngroup, rank, rank2, descr->segment_start, descr->segment_end, descr->segment_leader_rank);
    }
}

static void
_destroy_segment_group(struct SegmentGroupDescr * descr)
{

    MPI_Comm_free(&descr->Segment);
    MPI_Comm_free(&descr->Group);
    MPI_Comm_free(&descr->Leader);
}

static int
_aggregated(
            BigBlock * block,
            BigBlockPtr * ptr,
            ptrdiff_t offset, /* offset of the entire comm */
            size_t localsize,
            BigArray * array,
            int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array),
            int root,
            MPI_Comm comm);

static int
_throttle_action(MPI_Comm comm, int concurrency, BigBlock * block,
    BigBlockPtr * ptr,
    BigArray * array,
    int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array)
)
{
    int i;
    int ThisTask, NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    struct SegmentGroupDescr seggrp[1];

    if(concurrency <= 0) {
        concurrency = NTask;
    }

    size_t avgsegsize;
    size_t localsize = array->dims[0];
    size_t myoffset;
    size_t * sizes = malloc(sizeof(sizes[0]) * NTask);

    size_t totalsize = _collect_sizes(localsize, sizes, &myoffset, comm);

    /* try to create as many segments as number of groups (thus one segment per group) */
    avgsegsize = totalsize / concurrency;

    if(avgsegsize <= 0) avgsegsize = 1;

    /* no segment shall exceed the memory bound set by maxsegsize, since it will be collected to a single rank */
    if(avgsegsize > _BigFileAggThreshold) avgsegsize = _BigFileAggThreshold;

    _create_segment_group(seggrp, sizes, avgsegsize, concurrency, comm);

    free(sizes);

    int rt = 0;
    int segment;

    for(segment = seggrp->segment_start;
        segment < seggrp->segment_end;
        segment ++) {

        MPI_Barrier(seggrp->Group);

        if(0 != (rt = big_file_mpi_broadcast_anyerror(rt, seggrp->Group))) {
            /* failed , abort. */
            continue;
        } 
        if(seggrp->ThisSegment != segment) continue;

        /* use the offset on the first task in the SegGroup */
        size_t offset = myoffset;
        MPI_Bcast(&offset, 1, MPI_LONG, 0, seggrp->Segment);

        rt = _aggregated(block, ptr, offset, localsize, array, action, seggrp->segment_leader_rank, seggrp->Segment);

    }

    if(0 == (rt = big_file_mpi_broadcast_anyerror(rt, comm))) {
        /* no errors*/
        big_block_seek_rel(block, ptr, totalsize);
    }

    _destroy_segment_group(seggrp);
    return rt;
}

static int
_aggregated(
            BigBlock * block,
            BigBlockPtr * ptr,
            ptrdiff_t offset, /* offset of the entire comm */
            size_t localsize, /* offset of the entire comm */
            BigArray * array,
            int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array),
            int root,
            MPI_Comm comm)
{
    size_t elsize = big_file_dtype_itemsize(block->dtype) * block->nmemb;

    /* This will aggregate to the root and write */
    BigBlockPtr ptr1[1];
    /* use memcpy because older compilers doesn't like *ptr assignments */
    memcpy(ptr1, ptr, sizeof(BigBlockPtr));

    int i;
    int e = 0;
    int rank;
    int nrank;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nrank);

    BigArray garray[1], larray[1];
    BigArrayIter iarray[1], ilarray[1];
    void * lbuf = malloc(elsize * localsize);
    void * gbuf = NULL;

    int recvcounts[nrank];
    int recvdispls[nrank + 1];

    recvdispls[0] = 0;
    recvcounts[rank] = localsize;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvcounts, 1, MPI_INT, comm);

    int grouptotalsize = localsize;

    MPI_Allreduce(MPI_IN_PLACE, &grouptotalsize, 1, MPI_INT, MPI_SUM, comm);

    for(i = 0; i < nrank; i ++) {
        recvdispls[i + 1] = recvdispls[i] + recvcounts[i];
    }

    MPI_Datatype mpidtype;
    MPI_Type_contiguous(elsize, MPI_BYTE, &mpidtype);
    MPI_Type_commit(&mpidtype);

    big_array_init(larray, lbuf, block->dtype, 2, (size_t[]){localsize, block->nmemb}, NULL);

    big_array_iter_init(iarray, array);
    big_array_iter_init(ilarray, larray);

    if(rank == root) {
        gbuf = malloc(grouptotalsize * elsize);
        big_array_init(garray, gbuf, block->dtype, 2, (size_t[]){grouptotalsize, block->nmemb}, NULL);
    }

    if(action == big_block_write) {
        _dtype_convert(ilarray, iarray, localsize * block->nmemb);
        MPI_Gatherv(lbuf, recvcounts[rank], mpidtype,
                    gbuf, recvcounts, recvdispls, mpidtype, root, comm);
    }
    if(rank == root) {
        big_block_seek_rel(block, ptr1, offset);
        e = action(block, ptr1, garray);
    }
    if(action == big_block_read) {
        MPI_Scatterv(gbuf, recvcounts, recvdispls, mpidtype,
                    lbuf, localsize, mpidtype, root, comm);
        _dtype_convert(iarray, ilarray, localsize * block->nmemb);
    }

    if(rank == root) {
        free(gbuf);
    }
    free(lbuf);

    MPI_Type_free(&mpidtype);

    return big_file_mpi_broadcast_anyerror(e, comm);
}

int
big_block_mpi_write(BigBlock * block, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    int rt = _throttle_action(comm, concurrency, block, ptr, array, big_block_write);
    return rt;
}

int
big_block_mpi_read(BigBlock * block, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    int rt = _throttle_action(comm, concurrency, block, ptr, array, big_block_read);
    return rt;
}
