#include <string.h>
#include "pmpfft.h"
#include "msg.h"
#include <signal.h>

static void get_position(void * pdata, ptrdiff_t index, double pos[3]) {
    PMStore * p = (PMStore *)pdata;
    pos[0] = p->x[index][0];
    pos[1] = p->x[index][1];
    pos[2] = p->x[index][2];
}
static size_t pack(void * pdata, ptrdiff_t index, void * buf, int flags) {
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(ptr) memcpy(&ptr[s], &p->field[index], sizeof(p->field[0])); \
        s += sizeof(p->field[0]); \
        flags &= ~f; \
    }
    #define DISPATCHC(f, field, c) \
    if(HAS(flags, f)) { \
        if(ptr) memcpy(&ptr[s], &p->field[index][c], sizeof(p->field[0][0])); \
        s += sizeof(p->field[0][0]); \
        flags &= ~f; \
    }

    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)

    /* components */
    DISPATCHC(PACK_ACC_X, acc, 0)
    DISPATCHC(PACK_ACC_Y, acc, 1)
    DISPATCHC(PACK_ACC_Z, acc, 2)
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)

    #undef DISPATCH
    #undef DISPATCHC
    if(flags != 0) {
        msg_abort(-1, "Runtime Error, unknown packing field.\n");
    }
    return s;
}
static void unpack(void * pdata, ptrdiff_t index, void * buf, int flags) {
    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;

    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        memcpy(&p->field[index], &ptr[s], sizeof(p->field[0])); \
        s += sizeof(p->field[0]); \
        flags &= ~f; \
    }
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    #undef DISPATCH
    if(flags != 0) {
        msg_abort(-1, "Runtime Error, unknown unpacking field.\n");
    }
}
static void reduce(void * pdata, ptrdiff_t index, void * buf, int flags) {
    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;

    #define DISPATCHC(f, field, i) \
        if(HAS(flags, f)) { \
            p->field[index][i] += * ((typeof(p->field[index][i])*) &ptr[s]); \
            s += sizeof(p->field[index][i]); \
            flags &= ~f; \
        }
    DISPATCHC(PACK_ACC_X, acc, 0);
    DISPATCHC(PACK_ACC_Y, acc, 1);
    DISPATCHC(PACK_ACC_Z, acc, 2);
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)
    #undef DISPATCHC
    if(flags != 0) {
        msg_abort(-1, "Runtime Error, unknown unpacking field.\n");
    }
}
static struct {
    void * p;
    size_t s;
} AllocTable[128];
size_t used_bytes = 0;

int NAllocTable = 0;
static void * malloczero(size_t s) {
    used_bytes += s;
    void * p = pfft_malloc(s);
    AllocTable[NAllocTable].s = s;
    AllocTable[NAllocTable].p = p;
    NAllocTable ++;
    return p;
}

static void myfree(void * p) {
    static size_t max_used_bytes = 0;
    if(used_bytes > max_used_bytes) {
        max_used_bytes = used_bytes;
        msg_printf(info, "Peak memory usage on rank 0: %td bytes\n", max_used_bytes);
    }
    NAllocTable --;
    if(AllocTable[NAllocTable].p != p) {
        msg_abort(-1, "Allocation and Free is mismatched.\n");
    }
    size_t s = AllocTable[NAllocTable].s;
    used_bytes -= s;
    pfft_free(p);
}

void pm_store_alloc_bare(PMStore * p, size_t np_upper) {
    pm_store_init(p);
    p->x = p->iface.malloc(sizeof(p->x[0]) * np_upper);
    p->v = p->iface.malloc(sizeof(p->v[0]) * np_upper);
    p->id = p->iface.malloc(sizeof(p->id[0]) * np_upper);
    p->acc = NULL;
    p->dx1 = NULL;
    p->dx2 = NULL;
    p->np = 0; 
    p->np_upper = np_upper;
};

void pm_store_init(PMStore * p) {
    memset(p, 0, sizeof(p[0]));
    p->iface.malloc = malloczero;
    p->iface.free = myfree;
    p->iface.pack = pack;
    p->iface.unpack = unpack;
    p->iface.reduce = reduce;
    p->iface.get_position = get_position;
    p->iface.AllAttributes = PACK_ALL;
}
void pm_store_alloc(PMStore * p, size_t np_upper) {
    pm_store_alloc_bare(p, np_upper);

    p->acc = p->iface.malloc(sizeof(p->acc[0]) * np_upper);
    p->dx1 = p->iface.malloc(sizeof(p->dx1[0]) * np_upper);
    p->dx2 = p->iface.malloc(sizeof(p->dx2[0]) * np_upper);
};

size_t 
pm_store_alloc_evenly(PMStore * p, size_t np_total, double alloc_factor, MPI_Comm comm) 
{
    /* allocate for np_total cross all */
    int NTask;
    MPI_Comm_size(comm, &NTask);
    pm_store_alloc(p, (size_t)(1.0 * np_total / NTask * 2));
}

void pm_store_destroy(PMStore * p) {
    if(p->dx2) p->iface.free(p->dx2);
    if(p->dx1) p->iface.free(p->dx1);
    if(p->acc) p->iface.free(p->acc);
    p->iface.free(p->id);
    p->iface.free(p->v);
    p->iface.free(p->x);
    msg_printf(info, "Remainging allocation %td bytes, %d chunks\n", used_bytes, NAllocTable);
}

void pm_store_read(PMStore * p, char * datasource) {
    /* parse data soure and read */
}

void pm_store_write(PMStore * p, char * datasource) {
    /* parse data soure and write */
}

static void permute(void * data, int np, size_t elsize, int * ind) {
    void * tmp = malloc(elsize * np);
    int i;
    for(i = 0; i < np; i ++) {
        memcpy(((char*) tmp) + i * elsize, ((char*) data) + ind[i] * elsize, elsize);
    }
    memcpy(data, tmp, np * elsize);
    free(tmp);
}

static void pm_store_permute(PMStore * p, int * ind) {
    permute(p->x, p->np, sizeof(p->x[0]), ind);
    permute(p->v, p->np, sizeof(p->v[0]), ind);
    permute(p->id, p->np, sizeof(p->id[0]), ind);
    if(p->acc)
        permute(p->acc, p->np, sizeof(p->acc[0]), ind);
    if(p->dx1)
        permute(p->dx1, p->np, sizeof(p->dx1[0]), ind);
    if(p->dx2)
        permute(p->dx2, p->np, sizeof(p->dx2[0]), ind);
}

typedef int (pm_store_target_func)(void * pdata, ptrdiff_t index, void * data);

void pm_store_wrap(PMStore * p, double BoxSize[3]) {
    int i;
    int d;
    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            while(p->x[i][d] < 0) p->x[i][d] += BoxSize[d];
            while(p->x[i][d] >= BoxSize[d]) p->x[i][d] -= BoxSize[d];
        }
    } 
}

void pm_store_decompose(PMStore * p, pm_store_target_func target_func, void * data, MPI_Comm comm) {
    int * target = p->iface.malloc(sizeof(int) * p->np);
    int NTask, ThisTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    ptrdiff_t i;
    for(i = 0; i < p->np; i ++) {
        target[i] = target_func(p, i, data);
        if(ThisTask == target[i]) {
            target[i] = 0;
        } else {
            target[i] ++;
        }
    }
    /* do a bincount; offset by -1 because -1 is for self */
    int * count = calloc(NTask + 1, sizeof(int));
    int * offsets = calloc(NTask + 1, sizeof(int));
    int * sendcount = count + 1;
    int * recvcount = malloc(sizeof(int) * (NTask));
    int * recvoffset = malloc(sizeof(int) * (NTask));
    int * sendoffset = malloc(sizeof(int) * (NTask));

    for(i = 0; i < p->np; i ++) {
        count[target[i]] ++;
    }
    cumsum(offsets, count, NTask + 1);

    int * arg = p->iface.malloc(sizeof(int) * p->np);
    for(i = 0; i < p->np; i ++) {
        int offset = offsets[target[i]] ++;
        arg[offset] = i;
    }
    
    pm_store_permute(p, arg);

    p->iface.free(arg);
    p->iface.free(target);

    MPI_Alltoall(sendcount, 1, MPI_INT, 
                 recvcount, 1, MPI_INT, 
                 comm);

    size_t Nsend = cumsum(sendoffset, sendcount, NTask);
    size_t Nrecv = cumsum(recvoffset, recvcount, NTask);
    size_t elsize = p->iface.pack(NULL, 0, NULL, p->iface.AllAttributes);

    void * send_buffer = p->iface.malloc(elsize * Nsend);
    void * recv_buffer = p->iface.malloc(elsize * Nrecv);

    p->np -= Nsend;
    for(i = 0; i < Nsend; i ++) {
        p->iface.pack(p, i + p->np, (char*) send_buffer + i * elsize, p->iface.AllAttributes);
    }

    MPI_Datatype PTYPE;
    MPI_Type_contiguous(elsize, MPI_BYTE, &PTYPE);
    MPI_Type_commit(&PTYPE);
    MPI_Alltoallv_sparse(
            send_buffer, sendcount, sendoffset, PTYPE,
            recv_buffer, recvcount, recvoffset, PTYPE,
            comm);
    MPI_Type_free(&PTYPE);
    
    for(i = 0; i < Nrecv; i ++) {
        p->iface.unpack(p, i + p->np, (char*) recv_buffer + i * elsize, p->iface.AllAttributes);
    }

    p->np += Nrecv;

    p->iface.free(recv_buffer);
    p->iface.free(send_buffer);
    free(recvcount);
    free(recvoffset);
    free(sendoffset);
    free(offsets);
    free(count);
}

