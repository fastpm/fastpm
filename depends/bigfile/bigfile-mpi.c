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
_assign_colors(size_t glocalsize, size_t * sizes, int * color, int NTask)
{
    int i;

    size_t current_size = 0;
    int current_color = 0;
    for(i = 0; i < NTask; i ++) {
        current_size += sizes[i];
        color[i] = current_color;
        if(current_size > glocalsize) {
            current_size = 0;
            current_color ++;
        }
    }
    return color[NTask - 1] + 1;
}

static int
_aggregated(
            BigBlock * block,
            BigBlockPtr * ptr,
            ptrdiff_t offset, /* offset of the entire comm */
            size_t localsize,
            BigArray * array,
            int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array),
            MPI_Comm comm);

static int
_throttle_action(MPI_Comm comm, int concurrency, BigBlock * block,
    BigBlockPtr * ptr,
    BigArray * array,
    int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array)
)
{
    int ThisTask, NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    if(concurrency <= 0) {
        concurrency = NTask;
    }

    ptrdiff_t * sizes = malloc(sizeof(sizes[0]) * NTask);
    size_t * offsets = malloc(sizeof(offsets[0]) * (NTask + 1));

    ptrdiff_t totalsize;

    sizes[ThisTask] = array->dims[0];

    MPI_Datatype MPI_PTRDIFFT;

    if(sizeof(ptrdiff_t) == sizeof(long)) {
        MPI_PTRDIFFT = MPI_LONG;
    } else if(sizeof(ptrdiff_t) == sizeof(int)) {
            MPI_PTRDIFFT = MPI_INT;
    } else {
        /* Weird architecture indeed. */
        abort();
    }


    MPI_Allreduce(&sizes[ThisTask], &totalsize, 1, MPI_PTRDIFFT, MPI_SUM, comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sizes, 1, MPI_PTRDIFFT, comm);

    /* divide into groups of this size */

    ptrdiff_t glocalsize = totalsize / concurrency + 1;

    if (glocalsize > _BigFileAggThreshold) {
        glocalsize = _BigFileAggThreshold;
    }

    int * segid = malloc(sizeof(segid[0]) * NTask);
    int Nsegid = _assign_colors(glocalsize, sizes, segid, NTask);

    int * seg_writer = malloc(sizeof(seg_writer[0]) * Nsegid);
    int i;

    for(i = 0; i < Nsegid; i ++) {
        seg_writer[i] = i * concurrency / Nsegid;
    }

    int WriterID = seg_writer[segid[ThisTask]];
    MPI_Comm WriterGroup;
    int WriterGroupSize;
    int WriterGroupRank;

    MPI_Comm_split(comm, WriterID, ThisTask, &WriterGroup);

    MPI_Comm_size(WriterGroup, &WriterGroupSize);
    MPI_Comm_rank(WriterGroup, &WriterGroupRank);

    /* compute the unique segments on this writer */
    int Nseg_this_writer = 0;

    for(i = 0; i < Nsegid; i ++) {
        if(seg_writer[i] == WriterID) {
            Nseg_this_writer ++;
        }
    }

    int * unique_segs = malloc(sizeof(int) * Nseg_this_writer);

    Nseg_this_writer = 0;
    for(i = 0; i < Nsegid; i ++) {
        if(seg_writer[i] == WriterID) {
            unique_segs[Nseg_this_writer] = i;
            Nseg_this_writer ++;
        }
    }

    if (_big_file_mpi_verbose) {
        if (WriterGroupRank == 0) {
            printf("WriterGroup = %d WriterGroupSize = %d segs = [ ", WriterID, WriterGroupSize);
            for(i = 0; i < Nseg_this_writer; i ++) {
                printf("%d ", unique_segs[i]);
            }
            printf("] \n");
        }
    }

    /* offsets */
    offsets[0] = 0;
    for(i = 0; i < NTask; i ++) {
        offsets[i + 1] = offsets[i] + sizes[i];
    }

    MPI_Comm SegGroup;

    MPI_Comm_split(WriterGroup, segid[ThisTask], ThisTask, &SegGroup);

    int rt = 0;
    int active;
    for(active = 0; active < Nseg_this_writer; active++) {

        MPI_Barrier(WriterGroup);

        if(0 != (rt = big_file_mpi_broadcast_anyerror(rt, WriterGroup))) {
            /* failed , abort. */
            continue;
        } 
        if(segid[ThisTask] != unique_segs[active]) continue;

        /* offset on the root of SegGroup matters */
        rt = _aggregated(block, ptr, offsets[ThisTask], sizes[ThisTask], array, action, SegGroup);

    }

    free(segid);
    free(seg_writer);
    free(unique_segs);
    free(sizes);
    free(offsets);

    MPI_Comm_free(&SegGroup);
    MPI_Comm_free(&WriterGroup);

    if(0 == (rt = big_file_mpi_broadcast_anyerror(rt, comm))) {
        /* no errors*/
        big_block_seek_rel(block, ptr, totalsize);
    }
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

    if(rank == 0) {
        gbuf = malloc(grouptotalsize * elsize);
        big_array_init(garray, gbuf, block->dtype, 2, (size_t[]){grouptotalsize, block->nmemb}, NULL);
    }

    if(action == big_block_write) {
        _dtype_convert(ilarray, iarray, localsize * block->nmemb);
        MPI_Gatherv(lbuf, recvcounts[rank], mpidtype,
                    gbuf, recvcounts, recvdispls, mpidtype, 0, comm);
    }
    if(rank == 0) {
        big_block_seek_rel(block, ptr1, offset);
        e = action(block, ptr1, garray);
    }
    if(action == big_block_read) {
        MPI_Scatterv(gbuf, recvcounts, recvdispls, mpidtype,
                    lbuf, localsize, mpidtype, 0, comm);
        _dtype_convert(iarray, ilarray, localsize * block->nmemb);
    }

    if(rank == 0) {
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
