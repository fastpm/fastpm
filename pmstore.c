#include <string.h>
#include "pmpfft.h"
#include "permute.c"
#include "msort.c"

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
    }

    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_ACC_X, acc[0])
    DISPATCH(PACK_ACC_Y, acc[1])
    DISPATCH(PACK_ACC_Z, acc[2])

    #undef DISPATCH
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
    }
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_ACC_X, acc[0])
    DISPATCH(PACK_ACC_Y, acc[1])
    DISPATCH(PACK_ACC_Z, acc[2])
    #undef DISPATCH
}
static void reduce(void * pdata, ptrdiff_t index, void * buf, int flags) {
    PMStore * p = (PMStore *)pdata;
    size_t s = 0;
    char * ptr = (char*) buf;

    #define DISPATCH(f, field) \
        if(HAS(flags, f)) { \
            p->field[index] += * ((typeof(p->field[index])*) &ptr[s]); \
            s += sizeof(p->field[index]); \
        }
    DISPATCH(PACK_ACC_X, acc[0]);
    DISPATCH(PACK_ACC_Y, acc[1]);
    DISPATCH(PACK_ACC_Z, acc[2]);
    #undef DISPATCH
}

static void * malloczero(size_t s) {
    return calloc(s, 1);
}
void pm_store_init_bare(PMStore * p, size_t np_upper) {
    p->iface.malloc = malloczero;
    p->iface.free = free;
    p->iface.pack = pack;
    p->iface.unpack = unpack;
    p->iface.reduce = reduce;
    p->iface.get_position = get_position;
    p->x = p->iface.malloc(sizeof(p->x[0]) * np_upper);
    p->v = p->iface.malloc(sizeof(p->v[0]) * np_upper);
    p->id = p->iface.malloc(sizeof(p->id[0]) * np_upper);
    p->acc[0] = NULL;
    p->acc[1] = NULL;
    p->acc[2] = NULL;
    p->dx1 = NULL;
    p->dx2 = NULL;
    p->np = 0; 
};

void pm_store_init(PMStore * p, size_t np_upper) {
    pm_store_init_bare(p, np_upper);

    p->acc[0] = p->iface.malloc(sizeof(p->acc[0][0]) * np_upper);
    p->acc[1] = p->iface.malloc(sizeof(p->acc[1][0]) * np_upper);
    p->acc[2] = p->iface.malloc(sizeof(p->acc[2][0]) * np_upper);
    p->dx1 = p->iface.malloc(sizeof(p->dx1[0]) * np_upper);
    p->dx2 = p->iface.malloc(sizeof(p->dx2[0]) * np_upper);
};

void pm_store_destroy(PMStore * p) {
    p->iface.free(p->x);
    p->iface.free(p->v);
    p->iface.free(p->id);
    if(p->acc[0]) p->iface.free(p->acc[0]);
    if(p->acc[1]) p->iface.free(p->acc[1]);
    if(p->acc[2]) p->iface.free(p->acc[2]);
    if(p->dx1) p->iface.free(p->dx1);
    if(p->dx2) p->iface.free(p->dx2);
}

void pm_store_read(PMStore * p, char * datasource) {
    /* parse data soure and read */
}

void pm_store_write(PMStore * p, char * datasource) {
    /* parse data soure and write */
}

void pm_store_permute(PMStore * p, int * ind) {
    permute(p->x, p->np, sizeof(p->x[0]), ind);
    permute(p->v, p->np, sizeof(p->v[0]), ind);
    permute(p->id, p->np, sizeof(p->id[0]), ind);
    if(p->acc[0])
        permute(p->acc[0], p->np, sizeof(p->acc[0]), ind);
    if(p->acc[1])
        permute(p->acc[1], p->np, sizeof(p->acc[1]), ind);
    if(p->acc[2])
        permute(p->acc[2], p->np, sizeof(p->acc[2]), ind);

    if(p->dx1)
        permute(p->dx1, p->np, sizeof(p->dx1[0]) * 3, ind);
    if(p->dx2)
        permute(p->dx2, p->np, sizeof(p->dx2[0]) * 3, ind);
}
