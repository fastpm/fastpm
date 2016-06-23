#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bigfile-mpi.h"
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
    int rt;

    if(comm == MPI_COMM_NULL) return 0;

    MPI_Comm_rank(comm, &rank);
    if(rank == 0) { 
        rt = big_block_create(bb, basename, dtype, nmemb, Nfile, fsize);
    }
    MPI_Bcast(&rt, 1, MPI_INT, 0, comm);
    if(rt) {
        big_file_mpi_broadcast_error(0, comm);
        return rt;
    }
    big_block_mpi_broadcast(bb, 0, comm);
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

typedef struct {
    MPI_Comm group;
    int GroupSize;
    int GroupRank;
    size_t offset;
    size_t totalsize;
} ThrottlePlan;

static int _throttle_plan_create(ThrottlePlan * plan, MPI_Comm comm, int concurrency, size_t localsize)
{
    int ThisTask;
    int NTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    int color = ThisTask * concurrency / NTask;
    MPI_Comm_split(MPI_COMM_WORLD, color, ThisTask, &plan->group);
    MPI_Comm_size(plan->group, &plan->GroupSize);
    MPI_Comm_rank(plan->group, &plan->GroupRank);
    MPI_Allreduce(MPI_IN_PLACE, &plan->GroupSize, 1, MPI_INT, MPI_MAX, comm);

    ptrdiff_t offsets[NTask + 1];
    ptrdiff_t sizes[NTask];
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

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sizes, 1, MPI_PTRDIFFT, comm);
    ptrdiff_t size = 0;
    offsets[0] = 0;
    int i;
    for(i = 0; i < NTask; i ++) {
        offsets[i + 1] = offsets[i] + sizes[i];
    }
    plan->offset = offsets[ThisTask];
    plan->totalsize = sizes[NTask];
    return 0;
}

static int _throttle_plan_destroy(ThrottlePlan * plan)
{
    MPI_Comm_free(&plan->group);
}

int big_block_mpi_write(BigBlock * bb, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    /* FIXME: add exceptions */
    ThrottlePlan plan;
    _throttle_plan_create(&plan, comm, concurrency, array->dims[0]);

    BigBlockPtr ptr1 = * ptr;

    /* TODO: aggregrate if the array is sufficiently small */
    int i;
    for(i = 0; i < plan.GroupSize; i ++) {
        MPI_Barrier(plan.group);
        if (i != plan.GroupRank) continue;
        big_block_seek_rel(bb, &ptr1, plan.offset);
        big_block_write(bb, &ptr1, array);
    }

    big_block_seek_rel(bb, ptr, plan.totalsize);

    _throttle_plan_destroy(&plan);
    return 0;
}

int big_block_mpi_read(BigBlock * bb, BigBlockPtr * ptr, BigArray * array, int concurrency, MPI_Comm comm)
{
    /* FIXME: add exceptions */
    ThrottlePlan plan;
    _throttle_plan_create(&plan, comm, concurrency, array->dims[0]);

    BigBlockPtr ptr1 = * ptr;
    /* TODO: aggregrate if the array is sufficiently small */
    int i;
    for(i = 0; i < plan.GroupSize; i ++) {
        MPI_Barrier(plan.group);
        if (i != plan.GroupRank) continue;
        big_block_seek_rel(bb, &ptr1, plan.offset);
        big_block_read(bb, &ptr1, array);
    }

    big_block_seek_rel(bb, ptr, plan.totalsize);

    _throttle_plan_destroy(&plan);
    return 0;
}
