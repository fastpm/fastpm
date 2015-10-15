#include <string.h>
#include "pmpfft.h"
#include "msg.h"
#include <signal.h>
#include "walltime.h"

static void get_position(void * pdata, ptrdiff_t index, double pos[3]) {
    PMStore * p = (PMStore *)pdata;
    pos[0] = p->x[index][0];
    pos[1] = p->x[index][1];
    pos[2] = p->x[index][2];
}

static size_t pack(void * pdata, ptrdiff_t index, void * buf, int flags) {
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            if(ptr) memcpy(&ptr[s], &p->field[index], sizeof(p->field[0])); \
            s += sizeof(p->field[0]); \
            flags &= ~f; \
        } \
    }
    #define DISPATCHC(f, field, c) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            if(ptr) memcpy(&ptr[s], &p->field[index][c], sizeof(p->field[0][0])); \
            s += sizeof(p->field[0][0]); \
            flags &= ~f; \
        } \
    }

    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_Q, q)
    DISPATCH(PACK_ACC, acc)

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
        if(p->field) { \
            memcpy(&p->field[index], &ptr[s], sizeof(p->field[0])); \
            s += sizeof(p->field[0]); \
            flags &= ~f; \
        } \
    }
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_Q, q)
    DISPATCH(PACK_ACC, acc)
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
        if(p->field) { \
            p->field[index][i] += * ((typeof(p->field[index][i])*) &ptr[s]); \
            s += sizeof(p->field[index][i]); \
            flags &= ~f; \
        } \
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

void pm_store_init(PMStore * p) {
    memset(p, 0, sizeof(p[0]));
    p->iface.malloc = malloczero;
    p->iface.free = myfree;
    p->iface.pack = pack;
    p->iface.unpack = unpack;
    p->iface.reduce = reduce;
    p->iface.get_position = get_position;
}
void pm_store_alloc(PMStore * p, size_t np_upper, int attributes) {
    pm_store_init(p);

    p->np = 0; 
    p->np_upper = np_upper;
    p->attributes = attributes;
    if(attributes & PACK_POS)
        p->x = p->iface.malloc(sizeof(p->x[0]) * np_upper);
    else
        p->x = NULL;

    if(attributes & PACK_Q)
        p->q = p->iface.malloc(sizeof(p->q[0]) * np_upper);
    else
        p->q = NULL;

    if(attributes & PACK_VEL)
        p->v = p->iface.malloc(sizeof(p->v[0]) * np_upper);
    else
        p->v = NULL;

    if(attributes & PACK_ID)
        p->id = p->iface.malloc(sizeof(p->id[0]) * np_upper);
    else
        p->id = NULL;

    if(attributes & PACK_ACC)
        p->acc = p->iface.malloc(sizeof(p->acc[0]) * np_upper);
    else
        p->acc = NULL;

    if(attributes & PACK_DX1)
        p->dx1 = p->iface.malloc(sizeof(p->dx1[0]) * np_upper);
    else
        p->dx1 = NULL;

    if(attributes & PACK_DX2)
        p->dx2 = p->iface.malloc(sizeof(p->dx2[0]) * np_upper);
    else
        p->dx2 = NULL;
};

size_t 
pm_store_alloc_evenly(PMStore * p, size_t np_total, int attributes, double alloc_factor, MPI_Comm comm) 
{
    /* allocate for np_total cross all */
    int NTask;
    MPI_Comm_size(comm, &NTask);
    pm_store_alloc(p, (size_t)(1.0 * np_total / NTask * 2), attributes);
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
    walltime_measure("/Decompose/Pre");
    pm_store_permute(p, arg);
    walltime_measure("/Decompose/Permute");

    MPI_Barrier(comm);
    walltime_measure("/Decompose/Wait");

    p->iface.free(arg);
    p->iface.free(target);

    MPI_Alltoall(sendcount, 1, MPI_INT, 
                 recvcount, 1, MPI_INT, 
                 comm);

    walltime_measure("/Decompose/All1");

    size_t Nsend = cumsum(sendoffset, sendcount, NTask);
    size_t Nrecv = cumsum(recvoffset, recvcount, NTask);
    size_t elsize = p->iface.pack(p, 0, NULL, p->attributes);

    void * send_buffer = p->iface.malloc(elsize * Nsend);
    void * recv_buffer = p->iface.malloc(elsize * Nrecv);

    p->np -= Nsend;
    for(i = 0; i < Nsend; i ++) {
        p->iface.pack(p, i + p->np, (char*) send_buffer + i * elsize, p->attributes);
    }


    walltime_measure("/Decompose/Pack");
    MPI_Datatype PTYPE;
    MPI_Type_contiguous(elsize, MPI_BYTE, &PTYPE);
    MPI_Type_commit(&PTYPE);
    MPI_Alltoallv_sparse(
            send_buffer, sendcount, sendoffset, PTYPE,
            recv_buffer, recvcount, recvoffset, PTYPE,
            comm);
    MPI_Type_free(&PTYPE);
    walltime_measure("/Decompose/All2");
    
    for(i = 0; i < Nrecv; i ++) {
        p->iface.unpack(p, i + p->np, (char*) recv_buffer + i * elsize, p->attributes);
    }
    walltime_measure("/Decompose/Unpack");

    p->np += Nrecv;

    p->iface.free(recv_buffer);
    p->iface.free(send_buffer);
    free(recvcount);
    free(recvoffset);
    free(sendoffset);
    free(offsets);
    free(count);
}

void 
pm_store_set_lagrangian_position(PMStore * p, PM * pm, double shift[3]) 
{
    /* fill pdata with a uniform grid, respecting domain given by pm */
    
    ptrdiff_t ind;
    int d;
    p->np = 1;
    for(d = 0; d < 3; d++) {
        p->np *= pm->IRegion.size[d];
    }
    if(p->np > p->np_upper) {
        msg_abort(-1, "Need %td particles; %td allocated\n", p->np, p->np_upper);
    }
    ptrdiff_t ptrmax = 0;
    ptrdiff_t i[3] = {0};
    ptrdiff_t ptr = 0;
    for(ind = 0; ind < pm->IRegion.total; ind ++, pm_inc_i_index(pm, i)){
        uint64_t id;
        /* avoid the padded region */
        if(i[2] >= pm->IRegion.size[2]) continue;

        id = 0;
        for(d = 0; d < 3; d ++) {
            id = id * pm->Nmesh[d] + i[d] + pm->IRegion.start[d];
        }
//        msg_aprintf(debug, "Creating particle at offset %td i = %td %td %td\n", ptr, i[0], i[1], i[2]);
        for(d = 0; d < 3; d ++) {
            p->x[ptr][d] = (i[d] + pm->IRegion.start[d]) * (pm->BoxSize[d] / pm->Nmesh[d]) + shift[d];

            p->id[ptr]  = id;

            /* set q if it is allocated. */
            if(p->q) p->q[ptr][d] = p->x[ptr][d];
        }
        if(ptr > ptrmax) ptrmax = ptr;
        ptr ++;
    }
    if(ptrmax + 1 != p->np) {
        msg_abort(-1, "This is an internal error, particle number mismatched with grid. %td != %td, allocsize=%td, shape=(%td %td %td)\n", 
            ptrmax + 1, p->np, pm->allocsize, 
            pm->IRegion.size[0],
            pm->IRegion.size[1],
            pm->IRegion.size[2]

            );
    }
}

