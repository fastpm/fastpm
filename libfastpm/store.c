#include <string.h>

#include <mpi.h>
#include <pfft.h>
#include <gsl/gsl_rng.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"

#define HAS(a, b) ((a & b) != 0)

void fastpm_store_get_position(FastPMStore * p, ptrdiff_t index, double pos[3])
{
    pos[0] = p->x[index][0];
    pos[1] = p->x[index][1];
    pos[2] = p->x[index][2];
}

void fastpm_store_get_lagrangian_position(FastPMStore * p, ptrdiff_t index, double pos[3])
{
    pos[0] = p->q[index][0];
    pos[1] = p->q[index][1];
    pos[2] = p->q[index][2];
}

static size_t pack(FastPMStore * p, ptrdiff_t index, void * buf, enum FastPMPackFields flags) {
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            if(ptr) { \
                VALGRIND_CHECK_MEM_IS_DEFINED(&p->field[index], sizeof(p->field[0])); \
                memcpy(&ptr[s], &p->field[index], sizeof(p->field[0])); \
            } \
            s += sizeof(p->field[0]); \
            flags &= ~f; \
        } \
    }
    #define DISPATCHC(f, field, c) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            if(ptr) { \
                VALGRIND_CHECK_MEM_IS_DEFINED(&p->field[index][c], sizeof(p->field[0][0])); \
                memcpy(&ptr[s], &p->field[index][c], sizeof(p->field[0][0])); \
            } \
            s += sizeof(p->field[0][0]); \
            flags &= ~f; \
        } \
    }

    size_t s = 0;
    char * ptr = (char*) buf;
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_MASK, mask)
    DISPATCH(PACK_LENGTH, length)
    DISPATCH(PACK_DENSITY, rho)
    DISPATCH(PACK_POTENTIAL, potential)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_Q, q)
    DISPATCH(PACK_AEMIT, aemit)
    DISPATCH(PACK_ACC, acc)
    DISPATCH(PACK_TIDAL, tidal)
    DISPATCH(PACK_FOF, fof)

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
    DISPATCHC(PACK_POS_X, x, 0)
    DISPATCHC(PACK_POS_Y, x, 1)
    DISPATCHC(PACK_POS_Z, x, 2)
    DISPATCHC(PACK_TIDAL_XX, tidal, 0)
    DISPATCHC(PACK_TIDAL_YY, tidal, 1)
    DISPATCHC(PACK_TIDAL_ZZ, tidal, 2)
    DISPATCHC(PACK_TIDAL_XY, tidal, 3)
    DISPATCHC(PACK_TIDAL_YZ, tidal, 4)
    DISPATCHC(PACK_TIDAL_ZX, tidal, 5)

    #undef DISPATCH
    #undef DISPATCHC
    if(flags != 0) {
        fastpm_raise(-1, "Runtime Error, unknown unpacking field: %08X\n", flags);
    }
    return s;
}
static void unpack(FastPMStore * p, ptrdiff_t index, void * buf, enum FastPMPackFields flags) {
    size_t s = 0;
    char * ptr = (char*) buf;

    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            VALGRIND_CHECK_MEM_IS_DEFINED(&ptr[s], sizeof(p->field[0])); \
            memcpy(&p->field[index], &ptr[s], sizeof(p->field[0])); \
            s += sizeof(p->field[0]); \
            flags &= ~f; \
        } \
    }
    #define DISPATCHC(f, field, c) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            VALGRIND_CHECK_MEM_IS_DEFINED(&ptr[s], sizeof(p->field[0][0])); \
            memcpy(&p->field[index][c], &ptr[s], sizeof(p->field[0][0])); \
            s += sizeof(p->field[0][0]); \
            flags &= ~f; \
        } \
    }
    DISPATCH(PACK_POS, x)
    DISPATCH(PACK_VEL, v)
    DISPATCH(PACK_ID, id)
    DISPATCH(PACK_MASK, mask)
    DISPATCH(PACK_LENGTH, length)
    DISPATCH(PACK_DENSITY, rho)
    DISPATCH(PACK_POTENTIAL, potential)
    DISPATCH(PACK_DX1, dx1)
    DISPATCH(PACK_DX2, dx2)
    DISPATCH(PACK_Q, q)
    DISPATCH(PACK_AEMIT, aemit)
    DISPATCH(PACK_ACC, acc)
    DISPATCH(PACK_TIDAL, tidal)
    DISPATCH(PACK_FOF, fof)

    DISPATCHC(PACK_ACC_X, acc, 0)
    DISPATCHC(PACK_ACC_Y, acc, 1)
    DISPATCHC(PACK_ACC_Z, acc, 2)
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)
    DISPATCHC(PACK_POS_X, x, 0)
    DISPATCHC(PACK_POS_Y, x, 1)
    DISPATCHC(PACK_POS_Z, x, 2)
    DISPATCHC(PACK_TIDAL_XX, tidal, 0)
    DISPATCHC(PACK_TIDAL_YY, tidal, 1)
    DISPATCHC(PACK_TIDAL_ZZ, tidal, 2)
    DISPATCHC(PACK_TIDAL_XY, tidal, 3)
    DISPATCHC(PACK_TIDAL_YZ, tidal, 4)
    DISPATCHC(PACK_TIDAL_ZX, tidal, 5)
    #undef DISPATCH
    #undef DISPATCHC
    if(flags != 0) {
        fastpm_raise(-1, "Runtime Error, unknown unpacking field: %08X\n", flags);
    }
}
static void reduce(FastPMStore * p, ptrdiff_t index, void * buf, enum FastPMPackFields flags) {
    size_t s = 0;
    char * ptr = (char*) buf;

    #define DISPATCHC(f, field, i) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            p->field[index][i] += * ((__typeof__(p->field[index][i])*) &ptr[s]); \
            s += sizeof(p->field[index][i]); \
            flags &= ~f; \
        } \
    }
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            p->field[index] += * ((__typeof__(p->field[index])*) &ptr[s]); \
            s += sizeof(p->field[index]); \
            flags &= ~f; \
        } \
    }
    DISPATCH(PACK_DENSITY, rho)
    DISPATCH(PACK_POTENTIAL, potential)
    DISPATCHC(PACK_ACC_X, acc, 0);
    DISPATCHC(PACK_ACC_Y, acc, 1);
    DISPATCHC(PACK_ACC_Z, acc, 2);
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)
    DISPATCHC(PACK_POS_X, x, 0)
    DISPATCHC(PACK_POS_Y, x, 1)
    DISPATCHC(PACK_POS_Z, x, 2)
    DISPATCHC(PACK_TIDAL_XX, tidal, 0)
    DISPATCHC(PACK_TIDAL_YY, tidal, 1)
    DISPATCHC(PACK_TIDAL_ZZ, tidal, 2)
    DISPATCHC(PACK_TIDAL_XY, tidal, 3)
    DISPATCHC(PACK_TIDAL_YZ, tidal, 4)
    DISPATCHC(PACK_TIDAL_ZX, tidal, 5)
    #undef DISPATCHC
    #undef DISPATCH
    if(flags != 0) {
        fastpm_raise(-1, "Runtime Error, unknown unpacking field. %d\n", flags);
    }
}
static double to_double(FastPMStore * p, ptrdiff_t index, enum FastPMPackFields flags) {
    double rt = 0.;
    #define DISPATCHC(f, field, i) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            rt = p->field[index][i]; \
            flags &= ~f; \
            goto byebye; \
        } \
    }
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            rt = p->field[index]; \
            flags &= ~f; \
            goto byebye; \
        } \
    }
    DISPATCH(PACK_DENSITY, rho)
    DISPATCH(PACK_POTENTIAL, potential)
    DISPATCHC(PACK_ACC_X, acc, 0);
    DISPATCHC(PACK_ACC_Y, acc, 1);
    DISPATCHC(PACK_ACC_Z, acc, 2);
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)
    DISPATCHC(PACK_POS_X, x, 0)
    DISPATCHC(PACK_POS_Y, x, 1)
    DISPATCHC(PACK_POS_Z, x, 2)
    DISPATCHC(PACK_TIDAL_XX, tidal, 0)
    DISPATCHC(PACK_TIDAL_YY, tidal, 1)
    DISPATCHC(PACK_TIDAL_ZZ, tidal, 2)
    DISPATCHC(PACK_TIDAL_XY, tidal, 3)
    DISPATCHC(PACK_TIDAL_YZ, tidal, 4)
    DISPATCHC(PACK_TIDAL_ZX, tidal, 5)
    #undef DISPATCH
    #undef DISPATCHC
byebye:
    if(flags != 0) {
        fastpm_raise(-1, "Runtime Error, unknown unpacking field. %d\n", flags);
        return 0;
    } else {
        return rt;
    }
}
static void from_double(FastPMStore * p, ptrdiff_t index, enum FastPMPackFields flags, double value) {
    #define DISPATCHC(f, field, i) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            p->field[index][i] = value; \
            flags &= ~f; \
            goto byebye; \
        } \
    }
    #define DISPATCH(f, field) \
    if(HAS(flags, f)) { \
        if(p->field) { \
            p->field[index] = value; \
            flags &= ~f; \
            goto byebye; \
        } \
    }
    DISPATCH(PACK_DENSITY, rho)
    DISPATCH(PACK_POTENTIAL, potential)
    DISPATCHC(PACK_ACC_X, acc, 0);
    DISPATCHC(PACK_ACC_Y, acc, 1);
    DISPATCHC(PACK_ACC_Z, acc, 2);
    DISPATCHC(PACK_DX1_X, dx1, 0)
    DISPATCHC(PACK_DX1_Y, dx1, 1)
    DISPATCHC(PACK_DX1_Z, dx1, 2)
    DISPATCHC(PACK_DX2_X, dx2, 0)
    DISPATCHC(PACK_DX2_Y, dx2, 1)
    DISPATCHC(PACK_DX2_Z, dx2, 2)
    DISPATCHC(PACK_POS_X, x, 0)
    DISPATCHC(PACK_POS_Y, x, 1)
    DISPATCHC(PACK_POS_Z, x, 2)
    DISPATCHC(PACK_TIDAL_XX, tidal, 0)
    DISPATCHC(PACK_TIDAL_YY, tidal, 1)
    DISPATCHC(PACK_TIDAL_ZZ, tidal, 2)
    DISPATCHC(PACK_TIDAL_XY, tidal, 3)
    DISPATCHC(PACK_TIDAL_YZ, tidal, 4)
    DISPATCHC(PACK_TIDAL_ZX, tidal, 5)
    #undef DISPATCH
    #undef DISPATCHC
byebye:
    if(flags != 0) {
        fastpm_raise(-1, "Runtime Error, unknown unpacking field. %d\n", flags);
    }
}

static ptrdiff_t
_alignsize(ptrdiff_t size)
{
    /* align sizes */
    return ((size + 1024) / 1024) * 1024;
}

void
fastpm_store_init_details(FastPMStore * p,
                  size_t np_upper,
                  enum FastPMPackFields attributes,
                  enum FastPMMemoryLocation loc,
                  const char * file,
                  const int line)
{
    memset(p, 0, sizeof(p[0]));
    p->mem = _libfastpm_get_gmem();
    p->pack = pack;
    p->unpack = unpack;
    p->reduce = reduce;
    p->get_position = fastpm_store_get_position;
    p->to_double = to_double;
    p->from_double = from_double;

    p->np = 0; 
    p->np_upper = np_upper;
    p->attributes = attributes;

    int it;
    p->base = NULL;

    #define SIZEIT(column, PACK) \
        if(it == 0) { \
            (size += ((attributes & PACK) != 0) * _alignsize(sizeof(p->column[0]) * np_upper)); \
        } else { \
            if(attributes & PACK) { \
                p->column = (void*) (((char*) p->base) + offset); \
                offset += _alignsize(sizeof(p->column[0]) * np_upper); \
            } else { \
                p->column = NULL; \
            } \
        }

    ptrdiff_t size = 0;
    ptrdiff_t offset = 0;
    for(it = 0; it < 2; it ++) {
        SIZEIT(q, PACK_Q);
        SIZEIT(x, PACK_POS);
        SIZEIT(v, PACK_VEL);
        SIZEIT(id, PACK_ID);
        SIZEIT(length, PACK_LENGTH);
        SIZEIT(fof , PACK_FOF);
        SIZEIT(acc, PACK_ACC);
        SIZEIT(dx1, PACK_DX1);
        SIZEIT(dx2, PACK_DX2);
        SIZEIT(aemit, PACK_AEMIT);
        SIZEIT(rho, PACK_DENSITY);
        SIZEIT(potential, PACK_POTENTIAL);
        SIZEIT(tidal, PACK_TIDAL);
        SIZEIT(mask, PACK_MASK);

        if(it == 0) {
            p->base = fastpm_memory_alloc_details(p->mem, "FastPMStore", size, loc, file, line);
        }
    };
}

size_t 
fastpm_store_init_evenly(FastPMStore * p, size_t np_total, enum FastPMPackFields attributes, double alloc_factor, MPI_Comm comm) 
{
    /* allocate for np_total cross all */
    int NTask;
    MPI_Comm_size(comm, &NTask);

    size_t np_upper = (size_t)(1.0 * np_total / NTask * alloc_factor);

    MPI_Bcast(&np_upper, 1, MPI_LONG, 0, comm);
    fastpm_store_init(p, np_upper, attributes, FASTPM_MEMORY_HEAP);
    return 0;
}

size_t
fastpm_store_get_np_total(FastPMStore * p, MPI_Comm comm)
{
    long long np = p->np;
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, comm);
    return np;
}

void 
fastpm_store_destroy(FastPMStore * p) 
{
    fastpm_memory_free(p->mem, p->base);
}

void fastpm_store_read(FastPMStore * p, char * datasource) {
    /* parse data soure and read */
}

void fastpm_store_write(FastPMStore * p, char * datasource) {
    /* parse data soure and write */
}

static void permute(void * data, int np, size_t elsize, int * ind) {
    void * tmp = malloc(elsize * np);
    if(!tmp) {
        fastpm_raise(-1, "No memory for permuting\n");
    }
    int i;
    for(i = 0; i < np; i ++) {
        memcpy(((char*) tmp) + i * elsize, ((char*) data) + ind[i] * elsize, elsize);
    }
    memcpy(data, tmp, np * elsize);
    free(tmp);
}

static void fastpm_store_permute(FastPMStore * p, int * ind) {
    if(p->x)
        permute(p->x, p->np, sizeof(p->x[0]), ind);
    if(p->v)
        permute(p->v, p->np, sizeof(p->v[0]), ind);
    if(p->id)
        permute(p->id, p->np, sizeof(p->id[0]), ind);
    if(p->mask)
        permute(p->mask, p->np, sizeof(p->mask[0]), ind);
    if(p->fof)
        permute(p->fof, p->np, sizeof(p->fof[0]), ind);
    if(p->potential)
        permute(p->potential, p->np, sizeof(p->potential[0]), ind);
    if(p->tidal)
        permute(p->tidal, p->np, sizeof(p->tidal[0]), ind);
    if(p->aemit)
        permute(p->aemit, p->np, sizeof(p->aemit[0]), ind);
    if(p->q)
        permute(p->q, p->np, sizeof(p->q[0]), ind);
    if(p->acc)
        permute(p->acc, p->np, sizeof(p->acc[0]), ind);
    if(p->dx1)
        permute(p->dx1, p->np, sizeof(p->dx1[0]), ind);
    if(p->dx2)
        permute(p->dx2, p->np, sizeof(p->dx2[0]), ind);
}

static int * __sort_by_id_arg;
static FastPMStore *  __sort_by_id_p;
int _sort_by_id_cmpfunc(const void * p1, const void * p2)
{
    const int * i1 = (const int*) p1;
    const int * i2 = (const int*) p2;
    
    int v1 = (__sort_by_id_p->id[*i1] < __sort_by_id_p->id[*i2]);
    int v2 = (__sort_by_id_p->id[*i1] > __sort_by_id_p->id[*i2]);

    return v2 - v1;
}

void fastpm_store_sort_by_id(FastPMStore * p)
{
    int * arg = fastpm_memory_alloc(p->mem, "Temp", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);
    int i;
    for(i = 0; i < p->np; i ++) {
        arg[i] = i;
    }
    /* FIXME: copy some version of qsort_r */
    __sort_by_id_arg = arg;
    __sort_by_id_p = p;
    qsort(arg, p->np, sizeof(arg[0]), _sort_by_id_cmpfunc);
    fastpm_store_permute(p, arg);
    fastpm_memory_free(p->mem, arg);
}

void 
fastpm_store_wrap(FastPMStore * p, double BoxSize[3])
{
    int i;
    int d;
    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            while(p->x[i][d] < 0) p->x[i][d] += BoxSize[d];
            while(p->x[i][d] >= BoxSize[d]) p->x[i][d] -= BoxSize[d];
        }
    } 
}

int
FastPMTargetPM (FastPMStore * p, ptrdiff_t i, PM * pm)
{
    double pos[3];
    p->get_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
}

/*
static void
checkx(FastPMStore * p)
{
    ptrdiff_t i;

    int j = 0;
    for(i = 0; i < p->np; i ++) {
        if(p->x[i][0] > 0) {
            j = 1;
        }
        if(p->x[i][1] > 0) {
            j = 1;
        }
        if(p->x[i][2] > 0) {
            j = 1;
        }
    }
    printf("%d\n", j);
}
*/

int
fastpm_store_decompose(FastPMStore * p,
    fastpm_store_target_func target_func,
    void * data, MPI_Comm comm)
{

    VALGRIND_CHECK_MEM_IS_DEFINED(p->x, sizeof(p->x[0]) * p->np);

    int * target = fastpm_memory_alloc(p->mem, "Target", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);
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

    int * arg = fastpm_memory_alloc(p->mem, "PermArg", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);
    for(i = 0; i < p->np; i ++) {
        int offset = offsets[target[i]] ++;
        arg[offset] = i;
    }

    fastpm_store_permute(p, arg);

    VALGRIND_CHECK_MEM_IS_DEFINED(p->x, sizeof(p->x[0]) * p->np);

    fastpm_memory_free(p->mem, arg);
    fastpm_memory_free(p->mem, target);

    MPI_Alltoall(sendcount, 1, MPI_INT, 
                 recvcount, 1, MPI_INT, 
                 comm);

    size_t Nsend = cumsum(sendoffset, sendcount, NTask);
    size_t Nrecv = cumsum(recvoffset, recvcount, NTask);
    size_t elsize = p->pack(p, 0, NULL, p->attributes);

    size_t max_size = p->np + Nrecv - Nsend;

    MPI_Allreduce(MPI_IN_PLACE, &max_size, 1, MPI_LONG, MPI_MAX, comm);

    if(max_size > p->np_upper) {
        fastpm_info("Need %td particles on most expensive rank; %td allocated\n", max_size, p->np_upper);
        goto fail_oom;
    }

    void * send_buffer = fastpm_memory_alloc(p->mem, "SendBuf", elsize * Nsend, FASTPM_MEMORY_HEAP);
    void * recv_buffer = fastpm_memory_alloc(p->mem, "RecvBuf", elsize * Nrecv, FASTPM_MEMORY_HEAP);

    p->np -= Nsend;

    for(i = 0; i < Nsend; i ++) {
        p->pack(p, i + p->np, (char*) send_buffer + i * elsize, p->attributes);
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
        p->unpack(p, i + p->np, (char*) recv_buffer + i * elsize, p->attributes);
    }

    p->np += Nrecv;

    fastpm_memory_free(p->mem, recv_buffer);
    fastpm_memory_free(p->mem, send_buffer);

    free(recvcount);
    free(recvoffset);
    free(sendoffset);
    free(offsets);
    free(count);
    return 0;

    fail_oom:
        free(recvcount);
        free(recvoffset);
        free(sendoffset);
        free(offsets);
        free(count);
    return -1;
}

void
fastpm_store_get_q_from_id(FastPMStore * p, uint64_t id, double q[3])
{
    ptrdiff_t pabs[3];

    int d;
    for(d = 0; d < 3; d++) {
        pabs[d] = id / p->q_strides[d];
        id -= pabs[d] * p->q_strides[d];
    }

    for(d = 0; d < 3; d++) {
        q[d] = pabs[d] * p->q_scale[d];
        q[d] += p->q_shift[d];
    }
}

void
fastpm_store_fill(FastPMStore * p, PM * pm, double * shift, ptrdiff_t * Nc)
{
    /* fill p with a uniform grid, respecting domain given by pm. use a subsample ratio. 
     * (every subsample grid points) */
    if(Nc == NULL) {
        Nc = pm_nmesh(pm);
    }
    int d;
    p->np = 1;
    for(d = 0; d < 3; d++) {
        int start = pm->IRegion.start[d] * Nc[d] / pm->Nmesh[d];
        int end = (pm->IRegion.start[d] + pm->IRegion.size[d]) * Nc[d] / pm->Nmesh[d];
        p->np *= end - start;
    }
    if(p->np > p->np_upper) {
        fastpm_raise(-1, "Need %td particles; %td allocated\n", p->np, p->np_upper);
    }
    ptrdiff_t ptr = 0;

    for(d = 0; d < 3; d ++) {
        if(shift) 
            p->q_shift[d] = shift[d];
        else
            p->q_shift[d] = 0;

        p->q_scale[d] = pm->BoxSize[d] / Nc[d];
    }

    p->q_strides[0] = Nc[1] * Nc[2];
    p->q_strides[1] = Nc[2];
    p->q_strides[2] = 1;

    PMXIter iter;
    for(pm_xiter_init(pm, &iter);
       !pm_xiter_stop(&iter);
        pm_xiter_next(&iter)){
        ptrdiff_t pabs_start[3];
        ptrdiff_t pabs_end[3];
        ptrdiff_t ii, jj, kk;

        for(d = 0; d < 3; d ++) {
            pabs_start[d] = iter.iabs[d] * Nc[d] / pm->Nmesh[d];
            pabs_end[d] = (iter.iabs[d] + 1) * Nc[d] / pm->Nmesh[d];
        }
        for(ii = pabs_start[0]; ii < pabs_end[0]; ii ++)
        for(jj = pabs_start[1]; jj < pabs_end[1]; jj ++)
        for(kk = pabs_start[2]; kk < pabs_end[2]; kk ++) {
            ptrdiff_t pabs[3] = {ii, jj, kk};

            uint64_t id = pabs[2] * p->q_strides[2] +
                         pabs[1] * p->q_strides[1] +
                         pabs[0] * p->q_strides[0] ;

            if(p->id) p->id[ptr] = id;
            if(p->mask) p->mask[ptr] = 0;

            fastpm_store_get_q_from_id(p, id, &p->x[ptr][0]);

            if(p->q) {
                /* set q if it is allocated. */
                for(d = 0; d < 3; d ++) {
                    p->q[ptr][d] = p->x[ptr][d];
                }
            }
            ptr ++;
        }
    }
    if(ptr != p->np) {
        fastpm_raise(-1, "This is an internal error, particle number mismatched with grid. %td != %td, allocsize=%td, shape=(%td %td %td)\n", 
            ptr, p->np, pm->allocsize,
            pm->IRegion.size[0],
            pm->IRegion.size[1],
            pm->IRegion.size[2]
            );
    }
    p->a_x = p->a_v = 0.;
}

void 
fastpm_store_summary(FastPMStore * p, double dx1[3], double dx2[3], MPI_Comm comm) 
{

    int d;
    for(d = 0; d < 3; d ++) {
        dx1[d] = 0;
        dx2[d] = 0;
    }
    ptrdiff_t i;

#pragma omp parallel for
    for(i = 0; i < p->np; i ++) {
        int d;
        for(d =0; d < 3; d++) {
#pragma omp atomic
            dx1[d] += p->dx1[i][d] * p->dx1[i][d];
#pragma omp atomic
            dx2[d] += p->dx2[i][d] * p->dx2[i][d];
        } 
    }
    uint64_t Ntot = p->np;

    MPI_Allreduce(MPI_IN_PLACE, dx1, 3, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, dx2, 3, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ntot,   1, MPI_LONG,  MPI_SUM, comm);
    for(d =0; d < 3; d++) {
        dx1[d] /= Ntot;
        dx1[d] = sqrt(dx1[d]);
        dx2[d] /= Ntot;
        dx2[d] = sqrt(dx2[d]);
    }

}

static void
_fastpm_store_copy(FastPMStore * p, ptrdiff_t start, FastPMStore * po, ptrdiff_t offset)
{
    if(p->np - start + offset > po->np_upper) {
        fastpm_raise(-1, "Not enough storage in target FastPMStore: asking for %td but has %td\n", p->np, po->np_upper);
    }

    ptrdiff_t ncopy = p->np - start;

    if(po->x) memcpy(&po->x[offset], &p->x[start], sizeof(p->x[0][0]) * 3 * ncopy);
    if(po->q) memcpy(&po->q[offset], &p->q[start], sizeof(p->q[0][0]) * 3 * ncopy);
    if(po->v) memcpy(&po->v[offset], &p->v[start], sizeof(p->v[0][0]) * 3 * ncopy);
    if(po->acc) memcpy(&po->acc[offset], &p->acc[start], sizeof(p->acc[0][0]) * 3 * ncopy);
    if(po->dx1) memcpy(&po->dx1[offset], &p->dx1[start], sizeof(p->dx1[0][0]) * 3 * ncopy);
    if(po->dx2) memcpy(&po->dx2[offset], &p->dx2[start], sizeof(p->dx2[0][0]) * 3 * ncopy);
    if(po->id) memcpy(&po->id[offset], &p->id[start], sizeof(p->id[0]) * ncopy);
    if(po->mask) memcpy(&po->mask[offset], &p->mask[start], sizeof(p->mask[0]) * ncopy);
    if(po->aemit) memcpy(&po->aemit[offset], &p->aemit[start], sizeof(p->aemit[0]) * ncopy);
    if(po->potential) memcpy(&po->potential[offset], &p->potential[start], sizeof(p->potential[0]) * ncopy);
    if(po->rho) memcpy(&po->rho[offset], &p->rho[start], sizeof(p->potential[0]) * ncopy);
    if(po->tidal) memcpy(&po->tidal[offset], &p->tidal[start], sizeof(p->tidal[0]) * ncopy);

    po->np = offset + ncopy;
    po->a_x = p->a_x;
    po->a_v = p->a_v;
    if(po != p) {
        memcpy(po->q_strides, p->q_strides, 3 * sizeof(p->q_strides[0]));
        memcpy(po->q_scale, p->q_scale, 3 * sizeof(p->q_scale[0]));
        memcpy(po->q_shift, p->q_shift, 3 * sizeof(p->q_shift[0]));
    }
}

void
fastpm_store_copy(FastPMStore * p, FastPMStore * po)
{
    _fastpm_store_copy(p, 0, po, 0);
}

void
fastpm_store_append(FastPMStore * p, FastPMStore * po)
{
    _fastpm_store_copy(p, 0, po, po->np);
}

void
fastpm_store_fill_subsample_mask(FastPMStore * p,
        double fraction,
        uint8_t * mask,
        MPI_Comm comm)
{
    gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    /* set uncorrelated seeds */
    double seed=1231584; //FIXME: set this properly.
    gsl_rng_set(random_generator, seed);
    int d;

    for(d = 0; d < ThisTask * 8; d++) {
        seed = 0x7fffffff * gsl_rng_uniform(random_generator);
    }

    gsl_rng_set(random_generator, seed);

    memset(mask, 0, p->np);

    ptrdiff_t i;
    for(i=0; i < p->np; i++) {
        double rand_i = gsl_rng_uniform(random_generator);
        int flag = fraction > 1 || rand_i <= fraction;
        mask[i] = flag;
    }

    gsl_rng_free(random_generator);
}

/*
 * Create a subsample, keeping only those with mask == True
 * */
void
fastpm_store_subsample(FastPMStore * p, uint8_t * mask, FastPMStore * po)
{
    ptrdiff_t i;
    ptrdiff_t j;
    j = 0;

    for(i = 0; i < p->np; i ++) {
        if(!mask[i]) continue;
        /* avoid memcpy of same address if we are doing subsample inplace */
        if(p != po || j != i) {
            if(po->x) memcpy(po->x[j], p->x[i], sizeof(p->x[0][0]) * 3);
            if(po->q) memcpy(po->q[j], p->q[i], sizeof(p->q[0][0]) * 3);
            if(po->v) memcpy(po->v[j], p->v[i], sizeof(p->v[0][0]) * 3);
            if(po->acc) memcpy(po->acc[j], p->acc[i], sizeof(p->acc[0][0]) * 3);
            if(po->dx1) memcpy(po->dx1[j], p->dx1[i], sizeof(p->dx1[0][0]) * 3);
            if(po->dx2) memcpy(po->dx2[j], p->dx2[i], sizeof(p->dx2[0][0]) * 3);
            if(po->id) po->id[j] = p->id[i];
            if(po->mask) po->mask[j] = p->mask[i];
            if(po->aemit) po->aemit[j] = p->aemit[i];
            if(po->potential) po->potential[j] = p->potential[i];
            if(po->tidal) memcpy(po->tidal[j], p->tidal[i], sizeof(p->tidal[0][0]) * 6);
        }
        j ++;
    }
    po->np = j;
    po->a_x = p->a_x;
    po->a_v = p->a_v;

    if(po != p) {
        memcpy(po->q_strides, p->q_strides, 3 * sizeof(p->q_strides[0]));
        memcpy(po->q_scale, p->q_scale, 3 * sizeof(p->q_scale[0]));
        memcpy(po->q_shift, p->q_shift, 3 * sizeof(p->q_shift[0]));
    }
}


static ptrdiff_t
binary_search(double foo, double a[], size_t n) {
    ptrdiff_t left = 0, right = n;
    ptrdiff_t mid;
    /* find the right side, because hist would count the number < edge, not <= edge */
    if(a[left] > foo) {
        return 0;
    }
    if(a[right - 1] <= foo) {
        return n;
    }
    while(right - left > 1) {
        mid = left + ((right - left - 1) >> 1);
        double pivot = a[mid];
        if(pivot > foo) {
            right = mid + 1;
            /* a[right - 1] > foo; */
        } else {
            left = mid + 1;
            /* a[left] <= foo; */
        }
    }
    return left;
}

void
fastpm_store_histogram_aemit(FastPMStore * store,
        ptrdiff_t * hist,
        double * edges,
        size_t nbins,
        MPI_Comm comm)
{
    ptrdiff_t i;

    memset(hist, 0, sizeof(ptrdiff_t) * (nbins + 2));

    for(i = 0; i < store->np; i ++) {
        int ibin = binary_search(store->aemit[i], edges, nbins + 1);
        hist[ibin] ++;
    }
    MPI_Allreduce(MPI_IN_PLACE, hist, nbins + 2, MPI_LONG, MPI_SUM, comm);
}
