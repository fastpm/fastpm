#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bigfile-mpi.h"

/* disable aggregation by default */
static size_t _BigFileAggThreshold = 0;

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
    int rt = big_file_open(bf, basename);
    MPI_Barrier(comm);
    return rt;
}

int big_file_mpi_create(BigFile * bf, const char * basename, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int rt = big_file_create(bf, basename);
    MPI_Barrier(comm);
    return rt;
}
int big_file_mpi_open_block(BigFile * bf, BigBlock * block, const char * blockname, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    char * basename = alloca(strlen(bf->basename) + strlen(blockname) + 128);
    sprintf(basename, "%s/%s/", bf->basename, blockname);
    return big_block_mpi_open(block, basename, comm);
}

int big_file_mpi_create_block(BigFile * bf, BigBlock * block, const char * blockname, const char * dtype, int nmemb, int Nfile, size_t size, MPI_Comm comm) {
    size_t fsize[Nfile];
    int i;
    for(i = 0; i < Nfile; i ++) {
        fsize[i] = size * (i + 1) / Nfile 
                 - size * (i) / Nfile;
    }
    if(comm == MPI_COMM_NULL) return 0;
    char * basename = alloca(strlen(bf->basename) + strlen(blockname) + 128);
    if(0 != _big_file_mksubdir_r(bf->basename, blockname)) return -1;
    sprintf(basename, "%s/%s/", bf->basename, blockname);
    return big_block_mpi_create(block, basename, dtype, nmemb, Nfile, fsize, comm);
}

int big_file_mpi_close(BigFile * bf, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rt = big_file_close(bf);
    MPI_Barrier(comm);
    return rt;
}

static int big_block_mpi_broadcast(BigBlock * bb, int root, MPI_Comm comm);
static void big_file_mpi_broadcast_error(int root, MPI_Comm comm);


int big_block_mpi_open(BigBlock * bb, const char * basename, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int rt;
    if(rank == 0) { 
        rt = big_block_open(bb, basename);
    }
    MPI_Bcast(&rt, 1, MPI_INT, 0, comm);
    if(rt) {
        big_file_mpi_broadcast_error(0, comm);
        return rt;
    }
    big_block_mpi_broadcast(bb, 0, comm);
    return 0;
}

int big_block_mpi_create(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, size_t fsize[], MPI_Comm comm) {
    int rank;
    int NTask;
    int rt;

    if(comm == MPI_COMM_NULL) return 0;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &rank);

    if(rank == 0) {
        rt = _big_block_create_internal(bb, basename, dtype, nmemb, Nfile, fsize);
    }
    MPI_Bcast(&rt, 1, MPI_INT, 0, comm);
    if(rt) {
        big_file_mpi_broadcast_error(0, comm);
        return rt;
    }
    big_block_mpi_broadcast(bb, 0, comm);

    int i;
    for(i = bb->Nfile * rank / NTask; i < bb->Nfile * (rank + 1) / NTask; i ++) {
        FILE * fp = _big_file_open_a_file(bb->basename, i, "w");
        if(fp == NULL) return -1;
        fclose(fp);
    }

    return 0;
}

int big_block_mpi_close(BigBlock * block, MPI_Comm comm) {
    if(comm == MPI_COMM_NULL) return 0;
    int rank;
    MPI_Comm_rank(comm, &rank);
    unsigned int * checksum = alloca(sizeof(int) * block->Nfile);
    MPI_Reduce(block->fchecksum, checksum, block->Nfile, MPI_UNSIGNED, MPI_SUM, 0, comm);
    int dirty;
    MPI_Reduce(&block->dirty, &dirty, 1, MPI_INT, MPI_LOR, 0, comm);
    int rt;
    if(rank == 0) {
        int i;
        big_block_set_dirty(block, dirty);
        for(i = 0; i < block->Nfile; i ++) {
            block->fchecksum[i] = checksum[i];
        }
    } else {
        /* only the root rank updates */
        big_block_set_dirty(block, 0);
        big_attrset_set_dirty(block->attrset, 0);
    }
    rt = big_block_close(block);
    if(rt) {
        return rt;
    }
    MPI_Barrier(comm);
    return 0;
}

static void big_file_mpi_broadcast_error(int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    char * error = big_file_get_error_message();
    int errorlen;
    if(rank == root) {
        errorlen = strlen(error);
    }
    MPI_Bcast(&errorlen, 1, MPI_INT, root, comm);
    if(rank != root) {
        error = alloca(errorlen + 1);
    }
    MPI_Bcast(error, errorlen + 1, MPI_BYTE, root, comm);
    if(rank != root) {
        big_file_set_error_message(error);
    }
}
static int big_block_mpi_broadcast(BigBlock * bb, int root, MPI_Comm comm) {
    ptrdiff_t i;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int lname = 0;
    void * attrpack;
    size_t attrpacksize = 0;
    if(rank == root) {
        lname = strlen(bb->basename);
        attrpack = big_attrset_pack(bb->attrset, &attrpacksize);
    }
    MPI_Bcast(&lname, 1, MPI_INT, root, comm);
    MPI_Bcast(bb, sizeof(BigBlock), MPI_BYTE, root, comm);
    MPI_Bcast(&attrpacksize, sizeof(size_t), MPI_BYTE, root, comm);
    if(rank != root) {
        bb->basename = calloc(lname + 1, 1);
        bb->fsize = calloc(bb->Nfile, sizeof(size_t));
        bb->foffset = calloc(bb->Nfile + 1, sizeof(size_t));
        bb->fchecksum = calloc(bb->Nfile, sizeof(int));
        attrpack = malloc(attrpacksize);
    }
    MPI_Bcast(attrpack, attrpacksize, MPI_BYTE, root, comm);
    if(rank != root) {
        bb->attrset = big_attrset_unpack(attrpack);
    }
    free(attrpack);
    MPI_Bcast(bb->basename, lname + 1, MPI_BYTE, root, comm);
    MPI_Bcast(bb->fsize, sizeof(ptrdiff_t) * bb->Nfile, MPI_BYTE, root, comm);
    MPI_Bcast(bb->fchecksum, sizeof(int) * bb->Nfile, MPI_BYTE, root, comm);
    MPI_Bcast(bb->foffset, sizeof(ptrdiff_t) * (bb->Nfile + 1), MPI_BYTE, root, comm);
    return 0;
}

typedef struct ThrottlePlan {
    MPI_Comm group;
    int GroupSize;
    int GroupRank;
    size_t offset;
    size_t totalsize;
    size_t localsize;
    size_t grouptotalsize;
    size_t elsize;
    MPI_Datatype mpidtype;
    BigBlock * block;
    int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array);
    int (*execute)(struct ThrottlePlan * plan, BigBlockPtr * ptr, BigArray * array);
} ThrottlePlan;

static int
_throttle_plan_execute_agg(ThrottlePlan * plan,
            BigBlockPtr * ptr, BigArray * array);
static int
_throttle_plan_execute_turns(ThrottlePlan * plan,
            BigBlockPtr * ptr, BigArray * array);
static int
_throttle_plan_create(ThrottlePlan * plan, MPI_Comm comm, int concurrency, BigBlock * block, size_t localsize,
    int (*action)(BigBlock * bb, BigBlockPtr * ptr, BigArray * array)
)
{
    int ThisTask, NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    if(concurrency <= 0) {
        concurrency = NTask;
    }

    int color = ThisTask * concurrency / NTask;
    MPI_Comm_split(MPI_COMM_WORLD, color, ThisTask, &plan->group);
    MPI_Comm_size(plan->group, &plan->GroupSize);
    MPI_Comm_rank(plan->group, &plan->GroupRank);
    MPI_Allreduce(MPI_IN_PLACE, &plan->GroupSize, 1, MPI_INT, MPI_MAX, comm);

    ptrdiff_t offsets[NTask + 1], sizes[NTask];
    ptrdiff_t grouptotalsize;

    MPI_Datatype MPI_PTRDIFFT;

    if(sizeof(ptrdiff_t) == sizeof(long)) {
        MPI_PTRDIFFT = MPI_LONG;
    } else if(sizeof(ptrdiff_t) == sizeof(int)) {
            MPI_PTRDIFFT = MPI_INT;
    } else {
        /* Weird architecture indeed. */
        abort();
    }
    sizes[ThisTask] = localsize;

    MPI_Allreduce(&localsize, &grouptotalsize, 1, MPI_PTRDIFFT, MPI_SUM, plan->group);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sizes, 1, MPI_PTRDIFFT, comm);
    offsets[0] = 0;
    int i;
    for(i = 0; i < NTask; i ++) {
        offsets[i + 1] = offsets[i] + sizes[i];
    }
    plan->offset = offsets[ThisTask];
    plan->totalsize = offsets[NTask];
    plan->localsize = localsize;
    plan->grouptotalsize = grouptotalsize;
    plan->elsize = dtype_itemsize(block->dtype) * block->nmemb;
    plan->action = action;
    plan->block = block;

    MPI_Type_contiguous(plan->elsize, MPI_BYTE, &plan->mpidtype);
    MPI_Type_commit(&plan->mpidtype);

    size_t gbytes = plan->grouptotalsize * plan->elsize;
    if(gbytes <= _BigFileAggThreshold) {
        plan->execute = _throttle_plan_execute_agg;
    } else {
        plan->execute = _throttle_plan_execute_turns;
    }
    return 0;
}

static int _throttle_plan_destroy(ThrottlePlan * plan)
{
    MPI_Comm_free(&plan->group);
    MPI_Type_free(&plan->mpidtype);
}

static int
_throttle_plan_execute_agg(ThrottlePlan * plan,
            BigBlockPtr * ptr, BigArray * array)
{
    /* This will aggregate to the root and write */
    BigBlockPtr ptr1[1] = {*ptr};
    int i;
    int e = 0;
    BigArray garray[1], larray[1];
    BigArrayIter iarray[1], ilarray[1];
    void * lbuf = malloc(plan->elsize * plan->localsize);
    void * gbuf = NULL;

    int recvcounts[plan->GroupSize];
    int recvdispls[plan->GroupSize + 1];

    recvdispls[0] = 0;
    recvcounts[plan->GroupRank] = plan->localsize;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvcounts, 1, MPI_INT, plan->group);

    for(i = 0; i < plan->GroupSize; i ++) {
        recvdispls[i + 1] = recvdispls[i] + recvcounts[i];
    }

    big_array_init(larray, lbuf, plan->block->dtype, 2, (size_t[]){plan->localsize, plan->block->nmemb}, NULL);

    big_array_iter_init(iarray, array);
    big_array_iter_init(ilarray, larray);

    if(plan->GroupRank == 0) {
        gbuf = malloc(plan->grouptotalsize * plan->elsize);
        big_array_init(garray, gbuf, plan->block->dtype, 2, (size_t[]){plan->grouptotalsize, plan->block->nmemb}, NULL);
    }

    if(plan->action == big_block_write) {
        _dtype_convert(ilarray, iarray, plan->localsize * plan->block->nmemb);
        MPI_Gatherv(lbuf, recvcounts[plan->GroupRank], plan->mpidtype,
                    gbuf, recvcounts, recvdispls, plan->mpidtype, 0, plan->group);
    }
    if(plan->GroupRank == 0) {
        big_block_seek_rel(plan->block, ptr1, plan->offset);
        e = plan->action(plan->block, ptr1, garray);
    }
    if(plan->action == big_block_read) {
        MPI_Scatterv(gbuf, recvcounts, recvdispls, plan->mpidtype,
                    lbuf, plan->localsize, plan->mpidtype, 0, plan->group);
        _dtype_convert(iarray, ilarray, plan->localsize * plan->block->nmemb);
    }
    if(plan->GroupRank == 0) {
        free(gbuf);
    }
    free(lbuf);

    big_block_seek_rel(plan->block, ptr, plan->totalsize);
    return e;
}

static int
_throttle_plan_execute_turns(ThrottlePlan * plan,
            BigBlockPtr * ptr, BigArray * array)
{
    int i;
    int e = 0;
    BigBlockPtr ptr1[1] = {*ptr};
    for(i = 0; i < plan->GroupSize; i ++) {
        MPI_Allreduce(MPI_IN_PLACE, &e, 1, MPI_INT, MPI_LOR, plan->group);
        if (e) continue;
        if (i != plan->GroupRank) continue;
        big_block_seek_rel(plan->block, ptr1, plan->offset);
        e = plan->action(plan->block, ptr1, array);
    }

    big_block_seek_rel(plan->block, ptr, plan->totalsize);
    return e;
}

int
big_block_mpi_write(BigBlock * block, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    /* FIXME: make the exception collective */
    ThrottlePlan plan[1];
    _throttle_plan_create(plan, comm, concurrency, block, array->dims[0], big_block_write);

    int rt = plan->execute(plan, ptr, array);

    _throttle_plan_destroy(plan);
    return rt;
}

int
big_block_mpi_read(BigBlock * block, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    /* FIXME: make the exception collective */
    ThrottlePlan plan[1];

    _throttle_plan_create(plan, comm, concurrency, block, array->dims[0], big_block_read);
    int rt = plan->execute(plan, ptr, array);

    _throttle_plan_destroy(plan);
    return rt;
}
