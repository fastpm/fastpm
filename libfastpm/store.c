#include <string.h>

#include <mpi.h>
#include <pfft.h>
#include <gsl/gsl_rng.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"

#define HAS(a, b) ((a & b) != 0)

static void
pack_any(FastPMStore * p, ptrdiff_t index, int ci, void * packed)
{
    size_t elsize = p->_column_info[ci].elsize;
    memcpy(packed, p->columns[ci] + index * elsize, elsize);
}

static void
unpack_any(FastPMStore * p, ptrdiff_t index, int ci, void * packed)
{
    size_t elsize = p->_column_info[ci].elsize;
    memcpy(p->columns[ci] + index * elsize, packed, elsize);
}

static void
pack_member_any(FastPMStore * p, ptrdiff_t index, int ci, int memb, void * packed)
{
    size_t nmemb = p->_column_info[ci].nmemb ;
    if(memb > nmemb) {
        fastpm_raise(-1, "memb %d greater than nmemb %d", memb, nmemb);
    }

    size_t elsize = p->_column_info[ci].elsize;
    size_t membsize = p->_column_info[ci].membsize;

    memcpy(packed, p->columns[ci] + index * elsize + membsize * memb, membsize);
}

static void
reduce_member_f4_add(FastPMStore * p, ptrdiff_t index, int ci, int memb, void * packed)
{
    size_t nmemb = p->_column_info[ci].nmemb ;
    if(memb > nmemb || memb < 0) {
        fastpm_raise(-1, "memb %d greater than nmemb %d", memb, nmemb);
    }

    size_t elsize = p->_column_info[ci].elsize;

    float * ptr = (float*) (p->columns[ci] + index * elsize);
    float * ptr_packed = (float*) packed;

    ptr[memb] += ptr_packed[0];
}

static double
to_double_f4 (FastPMStore * p, ptrdiff_t index, int ci, int memb)
{
    size_t nmemb = p->_column_info[ci].nmemb ;
    if(memb > nmemb) {
        fastpm_raise(-1, "memb %d greater than nmemb %d", memb, nmemb);
    }
    size_t elsize = p->_column_info[ci].elsize;

    float * ptr = (float*) (p->columns[ci] + index * elsize);

    return ptr[memb];
}

static void
from_double_f4 (FastPMStore * p, ptrdiff_t index, int ci, int memb, const double value)
{
    size_t nmemb = p->_column_info[ci].nmemb ;
    if(memb > nmemb) {
        fastpm_raise(-1, "memb %d greater than nmemb %d", memb, nmemb);
    }
    size_t elsize = p->_column_info[ci].elsize;

    float * ptr = (float*) (p->columns[ci] + index * elsize);

    ptr[memb] = value;
}


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

static ptrdiff_t
_alignsize(ptrdiff_t size)
{
    /* align sizes */
    return ((size + 1024) / 1024) * 1024;
}

void
fastpm_store_init_details(FastPMStore * p,
                  size_t np_upper,
                  FastPMColumnTags attributes,
                  enum FastPMMemoryLocation loc,
                  const char * file,
                  const int line)
{
    memset(p, 0, sizeof(p[0]));
    p->mem = _libfastpm_get_gmem();

    p->np = 0; 
    p->np_upper = np_upper;
    p->attributes = attributes;

    /* clear the column pointers. */
    memset(p->columns, 0, sizeof(p->columns));
    memset(p->_column_info, 0, sizeof(p->_column_info));

    int it;
    p->_base = NULL;

    /* first loop initialize the _column_info struct for corresponding items in the pointer union; 
     * second loop initialize the pointers */
    #define DEFINE_COLUMN(column, attr_, dtype_, nmemb_) \
        { \
            int ci = FASTPM_STORE_COLUMN_INDEX(column); \
            if(attr_ != (1 << ci)) { fastpm_raise(-1, "attr and column are out of order for %s\n", # column ); } \
            strcpy(p->_column_info[ci].dtype, dtype_);  \
            p->_column_info[ci].elsize = sizeof(p->column[0]);  \
            p->_column_info[ci].nmemb = nmemb_;  \
            p->_column_info[ci].membsize = sizeof(p->column[0]) / nmemb_;  \
            p->_column_info[ci].attribute = attr_;  \
            p->_column_info[ci].pack = pack_any;  \
            p->_column_info[ci].unpack = unpack_any;  \
            p->_column_info[ci].pack_member = pack_member_any; \
            p->_column_info[ci].reduce_member = NULL;  \
            p->_column_info[ci].from_double = NULL;  \
            p->_column_info[ci].to_double = NULL;  \
        }

    #define COLUMN_INFO(column) (p->_column_info[FASTPM_STORE_COLUMN_INDEX(column)])

    DEFINE_COLUMN(x, PACK_POS, "f8", 3);
    DEFINE_COLUMN(q, PACK_Q, "f4", 3);
    DEFINE_COLUMN(v, PACK_VEL, "f4", 3);
    DEFINE_COLUMN(acc, PACK_ACC, "f4", 3);
    DEFINE_COLUMN(dx1, PACK_DX1, "f4", 3);
    DEFINE_COLUMN(dx2, PACK_DX2, "f4", 3);
    DEFINE_COLUMN(aemit, PACK_AEMIT, "f4", 1);
    DEFINE_COLUMN(rho, PACK_DENSITY, "f4", 1);
    DEFINE_COLUMN(potential, PACK_POTENTIAL, "f4", 1);
    DEFINE_COLUMN(tidal, PACK_TIDAL, "f4", 6);
    DEFINE_COLUMN(id, PACK_ID, "i8", 6);
    DEFINE_COLUMN(mask, PACK_MASK, "i1", 1);
    DEFINE_COLUMN(minid, PACK_MINID, "i8", 1);
    DEFINE_COLUMN(task, PACK_TASK, "i4", 1);
    DEFINE_COLUMN(length, PACK_LENGTH, "i4", 1);
    DEFINE_COLUMN(rdisp, PACK_RDISP, "f4", 6);
    DEFINE_COLUMN(vdisp, PACK_VDISP, "f4", 6);
    DEFINE_COLUMN(rvdisp, PACK_RVDISP, "f4", 9);

    COLUMN_INFO(rho).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(rho).from_double = from_double_f4;
    COLUMN_INFO(rho).to_double = to_double_f4;
    COLUMN_INFO(acc).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(acc).from_double = from_double_f4;
    COLUMN_INFO(dx1).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(dx1).from_double = from_double_f4;
    COLUMN_INFO(dx2).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(dx2).from_double = from_double_f4;
    COLUMN_INFO(potential).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(potential).from_double = from_double_f4;
    COLUMN_INFO(tidal).reduce_member = reduce_member_f4_add;
    COLUMN_INFO(tidal).from_double = from_double_f4;

    ptrdiff_t size = 0;
    ptrdiff_t offset = 0;
    for(it = 0; it < 2; it ++) {
        int ci;
        for(ci = 0; ci < 32; ci ++) {
            if(it == 0) {
                size += ((attributes & p->_column_info[ci].attribute) != 0) * (_alignsize(p->_column_info[ci].elsize * np_upper)); \
            } else { \
                if(attributes & p->_column_info[ci].attribute) { \
                    p->columns[ci] = (void*) (((char*) p->_base) + offset); \
                    offset += _alignsize(p->_column_info[ci].elsize * np_upper); \
                } else { \
                    p->columns[ci] = NULL; \
                } \
            }


        }
        if(it == 0) {
            p->_base = fastpm_memory_alloc_details(p->mem, "FastPMStore", size, loc, file, line);
            /* zero out all memory */
            memset(p->_base, 0, size);
        }
    };
}

size_t 
fastpm_store_init_evenly(FastPMStore * p, size_t np_total, FastPMColumnTags attributes, double alloc_factor, MPI_Comm comm) 
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

size_t
fastpm_store_get_mask_sum(FastPMStore * p, MPI_Comm comm)
{
    long long np = 0;
    ptrdiff_t i;
    for(i = 0; i < p->np; i ++) {
        np += p->mask[i] != 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, comm);
    return np;
}

void
fastpm_packing_plan_init(FastPMPackingPlan * plan, FastPMStore * p, FastPMColumnTags attributes)
{
    int ci;
    int i = 0;
    plan->elsize = 0;
    plan->attributes = attributes;
    for (ci = 0; ci < 32; ci ++) {
        if (!(p->_column_info[ci].attribute & attributes)) continue;
        plan->_ci[i] = ci;
        plan->_offsets[ci] = plan->elsize;
        plan->elsize += p->_column_info[ci].elsize;
        plan->_column_info[ci] = p->_column_info[ci];
        i++;
    }
    plan->Ncolumns = i;
}

void
fastpm_packing_plan_pack(FastPMPackingPlan * plan,
            FastPMStore * p, ptrdiff_t i, void * packed)
{
    int t;
    for (t = 0; t < plan->Ncolumns; t ++) {
        int ci = plan->_ci[t];
        ptrdiff_t offset = plan->_offsets[ci];
        plan->_column_info[ci].pack(p, i, ci,
            ((char*) packed) + offset);
    }
}

void
fastpm_packing_plan_unpack(FastPMPackingPlan * plan,
            FastPMStore * p, ptrdiff_t i, void * packed)
{
    int t;
    for (t = 0; t < plan->Ncolumns; t ++) {
        int ci = plan->_ci[t];
        fastpm_packing_plan_unpack_ci(plan, ci, p, i, packed);
    }
}

/* unpack a single column from the offset in packed data. */
void
fastpm_packing_plan_unpack_ci(FastPMPackingPlan * plan, int ci,
            FastPMStore * p, ptrdiff_t i, void * packed)
{
    ptrdiff_t offset = plan->_offsets[ci];
    plan->_column_info[ci].unpack(p, i, ci,
        ((char*) packed) + offset);
}


int
fastpm_store_find_column_id(FastPMStore * p, FastPMColumnTags attribute)
{
    int ci;
    for (ci = 0; ci < 32; ci ++) {
        if (p->_column_info[ci].attribute == attribute) {
            return ci;
        }
    }
    fastpm_raise(-1, "Unknown column %x", attribute);
    return -1;
}

void 
fastpm_store_destroy(FastPMStore * p) 
{
    fastpm_memory_free(p->mem, p->_base);
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

void fastpm_store_permute(FastPMStore * p, int * ind)
{
    int c;
    for(c = 0; c < 32; c ++) {
        if(!p->columns[c]) continue;
        permute(p->columns[c], p->np, p->_column_info[c].elsize, ind);
    }
}


static FastPMStore *  _fastpm_store_sort_store;
static int (*_fastpm_store_sort_cmp_func)(const int i1, const int i2, FastPMStore * p);

int _sort_by_id_cmpfunc(const void * p1, const void * p2)
{
    const int * i1 = (const int*) p1;
    const int * i2 = (const int*) p2;

    return _fastpm_store_sort_cmp_func(*i1, *i2, _fastpm_store_sort_store);
}

int
FastPMLocalSortByID(const int i1,
                    const int i2,
                    FastPMStore * p)
{
    int v1 = (p->id[i1] < p->id[i2]);
    int v2 = (p->id[i1] > p->id[i2]);

    return v2 - v1;
}

/* sort a store locally with in the MPI rank. 
 *
 * cmp_func is called with three arguments.
 * */
void
fastpm_store_sort(FastPMStore * p,
        int (*cmp_func)(const int i1, const int i2, FastPMStore * p))
{
    int * arg = fastpm_memory_alloc(p->mem, "Temp", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);
    int i;
    for(i = 0; i < p->np; i ++) {
        arg[i] = i;
    }
    /* FIXME: copy some version of qsort_r */
    _fastpm_store_sort_store = p;
    _fastpm_store_sort_cmp_func = cmp_func;
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
    fastpm_store_get_position(p, i, pos);
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

    if(p->np > p->np_upper) {
        fastpm_raise(-1, "Particle buffer overrun detected np = %td > np_upper %td.\n", p->np, p->np_upper);
    }
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
    FastPMPackingPlan plan[1];

    fastpm_packing_plan_init(plan, p, p->attributes);

    size_t elsize = plan->elsize;

    volatile size_t neededsize = p->np + Nrecv - Nsend;

    if(neededsize > p->np_upper) {
        fastpm_ilog(INFO, "Need %td particles on rank %d; %td allocated\n", neededsize, ThisTask, p->np_upper);
    }

    if(MPIU_Any(comm, neededsize > p->np_upper)) {
        goto fail_oom;
    }

    void * send_buffer = fastpm_memory_alloc(p->mem, "SendBuf", elsize * Nsend, FASTPM_MEMORY_HEAP);
    void * recv_buffer = fastpm_memory_alloc(p->mem, "RecvBuf", elsize * Nrecv, FASTPM_MEMORY_HEAP);

    p->np -= Nsend;

    for(i = 0; i < Nsend; i ++) {
        fastpm_packing_plan_pack(plan, p, i + p->np, (char*) send_buffer + i * elsize);
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
        fastpm_packing_plan_unpack(plan, p, i + p->np, (char*) recv_buffer + i * elsize);
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
        pabs[d] = id / p->_q_strides[d];
        id -= pabs[d] * p->_q_strides[d];
    }

    for(d = 0; d < 3; d++) {
        q[d] = pabs[d] * p->_q_scale[d];
        q[d] += p->_q_shift[d];
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
            p->_q_shift[d] = shift[d];
        else
            p->_q_shift[d] = 0;

        p->_q_scale[d] = pm->BoxSize[d] / Nc[d];
    }

    p->_q_strides[0] = Nc[1] * Nc[2];
    p->_q_strides[1] = Nc[2];
    p->_q_strides[2] = 1;

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

            uint64_t id = pabs[2] * p->_q_strides[2] +
                         pabs[1] * p->_q_strides[1] +
                         pabs[0] * p->_q_strides[0] ;

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
_fastpm_store_copy(FastPMStore * p, ptrdiff_t start, FastPMStore * po, ptrdiff_t offset, size_t ncopy)
{
    if(ncopy + start > p->np) {
        fastpm_raise(-1, "Copy out of bounds from source FastPMStore: asking for %td but has %td\n", ncopy + start, p->np);
    }
    if(ncopy + offset > po->np_upper) {
        fastpm_raise(-1, "Not enough storage in target FastPMStore: asking for %td but has %td\n", p->np, po->np_upper);
    }

    int c;
    for(c = 0; c < 32; c ++) {
        if(!po->columns[c]) continue;
        size_t elsize = po->_column_info[c].elsize;

        memcpy(po->columns[c] + offset * elsize,
                p->columns[c] + start * elsize, elsize * ncopy);

    }
    po->np = offset + ncopy;
    po->a_x = p->a_x;
    po->a_v = p->a_v;
    if(po != p) {
        memcpy(po->_q_strides, p->_q_strides, 3 * sizeof(p->_q_strides[0]));
        memcpy(po->_q_scale, p->_q_scale, 3 * sizeof(p->_q_scale[0]));
        memcpy(po->_q_shift, p->_q_shift, 3 * sizeof(p->_q_shift[0]));
    }
}

void
fastpm_store_copy(FastPMStore * p, FastPMStore * po)
{
    _fastpm_store_copy(p, 0, po, 0, p->np);
}

void
fastpm_store_take(FastPMStore * p, ptrdiff_t i, FastPMStore * po, ptrdiff_t j)
{
    _fastpm_store_copy(p, i, po, j, 1);
}

void
fastpm_store_append(FastPMStore * p, FastPMStore * po)
{
    _fastpm_store_copy(p, 0, po, po->np, p->np);
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
        if(p == po && j == i) {j ++; continue; }

        int c;
        for(c = 0; c < 32; c ++) {
            if (!po->columns[c]) continue;

            size_t elsize = po->_column_info[c].elsize;
            memcpy(po->columns[c] + j * elsize, p->columns[c] + i * elsize, elsize);
        }

        j ++;
    }

    po->np = j;
    po->a_x = p->a_x;
    po->a_v = p->a_v;

    if(po != p) {
        memcpy(po->_q_strides, p->_q_strides, 3 * sizeof(p->_q_strides[0]));
        memcpy(po->_q_scale, p->_q_scale, 3 * sizeof(p->_q_scale[0]));
        memcpy(po->_q_shift, p->_q_shift, 3 * sizeof(p->_q_shift[0]));
    }
}


#if 0
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
#endif
/* this is cumulative */
void
fastpm_store_histogram_aemit_sorted(FastPMStore * store,
        int64_t * hist,
        double * aedges,
        size_t nedges,
        MPI_Comm comm)
{
    ptrdiff_t i;

    int64_t * hist1 = malloc(sizeof(hist1[0]) * (nedges + 1));

    memset(hist1, 0, sizeof(hist1[0]) * (nedges + 1));

#pragma omp parallel
    {
        /* FIXME: use standard reduction with OpenMP 4.7 */

        int64_t * hist2 = malloc(sizeof(hist2[0]) * (nedges + 1));
        memset(hist2, 0, sizeof(hist2[0]) * (nedges + 1));

        int iedge = 0;

        /* this works because openmp will send each thread at most 1
         * chunk with the static scheduling; thus a thread never sees
         * out of order aemit */
        #pragma omp for schedule(static)
        for(i = 0; i < store->np; i ++) {
            while(iedge < nedges && store->aemit[i] >= aedges[iedge]) iedge ++;

//           int ibin = binary_search(store->aemit[i], aedges, nedges);
//           if(ibin != iedge) abort();

            hist2[iedge] ++;
        }

        #pragma omp critical
        {
            for(i = 0; i < nedges + 1; i ++) {
                hist1[i] += hist2[i];
            }
        }

        free(hist2);
    }

    MPI_Allreduce(MPI_IN_PLACE, hist1, nedges + 1, MPI_LONG, MPI_SUM, comm);

    for(i = 0; i < nedges + 1; i ++) {
        hist[i] += hist1[i];
    }

    free(hist1);
}

