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

void
FastPMReduceOverwriteAny(FastPMStore * p, ptrdiff_t index, int ci, void * packed, void * userdata)
{
    size_t elsize = p->_column_info[ci].elsize;

    memcpy(p->columns[ci] + index * elsize, packed, elsize);
}

void
FastPMReduceAddFloat(FastPMStore * src, ptrdiff_t isrc, FastPMStore * dest, ptrdiff_t idest, int ci, void * userdata)
{
    size_t nmemb = dest->_column_info[ci].nmemb ;

    size_t elsize = dest->_column_info[ci].elsize;

    float * pdest = (float*) (dest->columns[ci] + idest * elsize);
    float * psrc  = (float*) (src ->columns[ci] + isrc * elsize);

    int d;
    for(d = 0; d < nmemb; d ++) {
        pdest[d] += psrc[d];
    }
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

static double
to_double_f8 (FastPMStore * p, ptrdiff_t index, int ci, int memb)
{
    size_t nmemb = p->_column_info[ci].nmemb ;
    if(memb > nmemb) {
        fastpm_raise(-1, "memb %d greater than nmemb %d", memb, nmemb);
    }
    size_t elsize = p->_column_info[ci].elsize;

    double * ptr = (double*) (p->columns[ci] + index * elsize);

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

const char *
fastpm_species_get_name(enum FastPMSpecies species)
{
    switch(species) {
        case FASTPM_SPECIES_BARYON:
            return "0";
        case FASTPM_SPECIES_CDM:
            return "1";
        case FASTPM_SPECIES_NCDM:
            return "2";
    }
    return "UNKNOWN";
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

double fastpm_store_get_mass(FastPMStore * p, ptrdiff_t index)
{
    /* total mass is the sum of the base and the extra */
    if(p->mass) {
        return p->meta.M0 + p->mass[index];
    } else {
        return p->meta.M0;
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
                  const char * name,
                  size_t np_upper,
                  FastPMColumnTags attributes,
                  enum FastPMMemoryLocation loc,
                  const char * file,
                  const int line)
{
    p->mem = _libfastpm_get_gmem();

    /* only set the name if name is not NULL; this is to allow setting the name before calling init.
     * */
    if(name) {
        strcpy(p->name, name);
    }

    p->attributes = attributes;

    p->np = 0;
    p->np_upper = np_upper;

    /* clear the column pointers. */
    memset(p->columns, 0, sizeof(p->columns));
    memset(p->_column_info, 0, sizeof(p->_column_info));

    /* clear the meta */
    memset(&p->meta, 0, sizeof(p->meta));
    int it;
    p->_base = NULL;

    /* first loop initialize the _column_info struct for corresponding items in the pointer union; 
     * second loop initialize the pointers */
    #define DEFINE_COLUMN(column, attr_, dtype_, nmemb_) \
        { \
            int ci = FASTPM_STORE_COLUMN_INDEX(column); \
            if(attr_ != (1 << ci)) { fastpm_raise(-1, "attr and column are out of order for %s\n", # column ); } \
            strcpy(p->_column_info[ci].dtype, dtype_);  \
            strcpy(p->_column_info[ci].name, # column); \
            p->_column_info[ci].elsize = sizeof(p->column[0]);  \
            p->_column_info[ci].nmemb = nmemb_;  \
            p->_column_info[ci].membsize = sizeof(p->column[0]) / nmemb_;  \
            p->_column_info[ci].attribute = attr_;  \
            p->_column_info[ci].pack = pack_any;  \
            p->_column_info[ci].unpack = unpack_any;  \
            p->_column_info[ci].from_double = NULL;  \
            p->_column_info[ci].to_double = NULL;  \
        }

    #define COLUMN_INFO(column) (p->_column_info[FASTPM_STORE_COLUMN_INDEX(column)])

    DEFINE_COLUMN(x, COLUMN_POS, "f8", 3);
    DEFINE_COLUMN(q, COLUMN_Q, "f4", 3);
    DEFINE_COLUMN(v, COLUMN_VEL, "f4", 3);
    DEFINE_COLUMN(acc, COLUMN_ACC, "f4", 3);
    DEFINE_COLUMN(dx1, COLUMN_DX1, "f4", 3);
    DEFINE_COLUMN(dx2, COLUMN_DX2, "f4", 3);
    DEFINE_COLUMN(dv1, COLUMN_DV1, "f4", 3);
    DEFINE_COLUMN(aemit, COLUMN_AEMIT, "f4", 1);
    DEFINE_COLUMN(rho, COLUMN_DENSITY, "f4", 1);
    DEFINE_COLUMN(potential, COLUMN_POTENTIAL, "f4", 1);
    DEFINE_COLUMN(tidal, COLUMN_TIDAL, "f4", 6);
    DEFINE_COLUMN(id, COLUMN_ID, "i8", 1);
    DEFINE_COLUMN(pgdc, COLUMN_PGDC, "f4", 3);
    DEFINE_COLUMN(mask, COLUMN_MASK, "i1", 1);
    DEFINE_COLUMN(minid, COLUMN_MINID, "i8", 1);
    DEFINE_COLUMN(task, COLUMN_TASK, "i4", 1);
    DEFINE_COLUMN(length, COLUMN_LENGTH, "i4", 1);
    DEFINE_COLUMN(rdisp, COLUMN_RDISP, "f4", 6);
    DEFINE_COLUMN(vdisp, COLUMN_VDISP, "f4", 6);
    DEFINE_COLUMN(rvdisp, COLUMN_RVDISP, "f4", 9);
    DEFINE_COLUMN(mass, COLUMN_MASS, "f4", 1);

    COLUMN_INFO(x).to_double = to_double_f8;
    COLUMN_INFO(v).to_double = to_double_f4;
    COLUMN_INFO(rho).to_double = to_double_f4;
    COLUMN_INFO(dx1).to_double = to_double_f4;
    COLUMN_INFO(dx2).to_double = to_double_f4;
    COLUMN_INFO(dv1).to_double = to_double_f4;
    COLUMN_INFO(acc).to_double = to_double_f4;
    COLUMN_INFO(mass).to_double = to_double_f4;

    COLUMN_INFO(rho).from_double = from_double_f4;
    COLUMN_INFO(acc).from_double = from_double_f4;
    COLUMN_INFO(pgdc).from_double = from_double_f4;
    COLUMN_INFO(dx1).from_double = from_double_f4;
    COLUMN_INFO(dx2).from_double = from_double_f4;
    COLUMN_INFO(dv1).from_double = from_double_f4;
    COLUMN_INFO(potential).from_double = from_double_f4;
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

void
fastpm_store_set_name(FastPMStore * p, const char * name)
{
    strcpy(p->name, name);
}

size_t
fastpm_store_init_evenly_details(FastPMStore * p,
    const char * name, size_t np_total, FastPMColumnTags attributes, double alloc_factor, MPI_Comm comm,
    const char * file,
    const int line)
{
    /* allocate for np_total cross all */
    /* name means name of species*/
    int NTask;
    MPI_Comm_size(comm, &NTask);

    size_t np_upper = (size_t)(1.0 * np_total / NTask * alloc_factor);

    MPI_Bcast(&np_upper, 1, MPI_LONG, 0, comm);
    fastpm_store_init_details(p, name, np_upper, attributes, FASTPM_MEMORY_HEAP, file, line);
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
        ptrdiff_t elsize = p->_column_info[ci].elsize;
        plan->elsize += elsize;
        plan->_column_info[ci] = p->_column_info[ci];
        i++;
    }
    /* Pad the elsize to 8 bytes. This ensures anything goes to the MPI wire is 8 byte aligned,
     * and minimizes the chances of hitting an implementation bug. */
    if (plan->elsize % 8 != 0) {
        plan->elsize += (8 - plan->elsize % 8);
    }
    plan->Ncolumns = i;
}

void
fastpm_packing_plan_pack(FastPMPackingPlan * plan,
            FastPMStore * p, ptrdiff_t i, void * packed)
{
    int t;
    memset(packed, 0, plan->elsize);
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
            double n = abs(p->x[i][d] / BoxSize[d]);

            double x1 = remainder(p->x[i][d], BoxSize[d]);

            while(x1 < 0) x1 += BoxSize[d];
            while(x1 > BoxSize[d]) x1 -= BoxSize[d];
            p->x[i][d] = x1;

            if(n > 10000) {
                double q[3];
                if(fastpm_store_has_q(p)) {
                    fastpm_store_get_q_from_id(p, p->id[i], q);
                }
                fastpm_raise(-1, "Particle at %g %g %g (q = %g %g %g) is too far from the bounds. Wrapping failed.\n", 
                        p->x[i][0],
                        p->x[i][1],
                        p->x[i][2],
                        q[0], q[1], q[2]
                );
            }
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

int
fastpm_store_decompose(FastPMStore * p,
    fastpm_store_target_func target_func,
    void * data, MPI_Comm comm)
{
    if(fastpm_store_get_np_total(p, comm) == 0) return 0 ;

    VALGRIND_CHECK_MEM_IS_DEFINED(p->x, sizeof(p->x[0]) * p->np);

    FastPMPackingPlan plan[1];

    fastpm_packing_plan_init(plan, p, p->attributes);

    size_t elsize = plan->elsize;

    size_t Nsend_limit = 1000 * 1024 * 1024 / elsize;     // this large number effectively prevents throttling 

    int NTask, ThisTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    if(p->np > p->np_upper) {
        fastpm_raise(-1, "Particle buffer overrun detected np = %td > np_upper %td.\n", p->np, p->np_upper);
    }

    /* do a bincount; offset by -1 because -1 is for self */
    int * count = calloc(NTask + 1, sizeof(int));
    int * offsets = calloc(NTask + 1, sizeof(int));
    int * sendcount = count + 1;
    int * recvcount = malloc(sizeof(int) * (NTask));
    int * recvoffset = malloc(sizeof(int) * (NTask));
    int * sendoffset = malloc(sizeof(int) * (NTask));

    int incomplete = 1;
    int iter = 0;
    /* terminate based on incomplete for throttling*/
    while(1) {
        incomplete = 0;
        int * target = fastpm_memory_alloc(p->mem, "Target", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);

        size_t Nsend_all = 0;
        ptrdiff_t i;
        for(i = 0; i < p->np; i ++) {
            target[i] = target_func(p, i, data);
            if(ThisTask == target[i]) {
                target[i] = -1;
            } else {
                Nsend_all++;
                /* Throttling: never send more than this many particles. */
                if(Nsend_all >= Nsend_limit) {
                    incomplete = 1;
                    target[i] = -1;
                } else {
                }
            }
        }

        for(i = 0; i < NTask + 1; i ++) {
            count[i] = 0;
            offsets[i] = 0;
        }
        for(i = 0; i < NTask; i ++) {
            sendoffset[i] = 0;
            recvoffset[i] = 0;
        }

        for(i = 0; i < p->np; i ++) {
            sendcount[target[i]] ++;
        }
        cumsum(offsets, count, NTask + 1);

        int * arg = fastpm_memory_alloc(p->mem, "PermArg", sizeof(int) * p->np, FASTPM_MEMORY_HEAP);
        for(i = 0; i < p->np; i ++) {
            int offset = offsets[target[i] + 1] ++;
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

        volatile size_t neededsize = p->np + Nrecv - Nsend;

        if(neededsize > p->np_upper) {
            fastpm_ilog(INFO, "Need %td particles on rank %d; %td allocated\n", neededsize, ThisTask, p->np_upper);
        }

        if(MPIU_Any(comm, neededsize > p->np_upper)) {
            goto fail_oom;
        }

        p->np -= Nsend;

        void * send_buffer = fastpm_memory_alloc(p->mem, "SendBuf", elsize * Nsend, FASTPM_MEMORY_HEAP);
        void * recv_buffer = fastpm_memory_alloc(p->mem, "RecvBuf", elsize * Nrecv, FASTPM_MEMORY_HEAP);

        {
            double nmin, nmax, nmean, nstd;

            MPIU_stats(comm, elsize * Nsend, "<>-s", &nmin, &nmax, &nmean, &nstd);
            fastpm_info("Send buffer size : min=%g max=%g mean=%g, std=%g bytes", nmin, nmax, nmean, nstd);
            MPIU_stats(comm, elsize * Nrecv, "<>-s", &nmin, &nmax, &nmean, &nstd);
            fastpm_info("Recv buffer size : min=%g max=%g mean=%g, std=%g bytes", nmin, nmax, nmean, nstd);
        }

        int ipar = 0;
        int j;
        for(i = 0; i < NTask; i ++) {
            for(j = sendoffset[i]; j < sendoffset[i] + sendcount[i]; j ++, ipar++) {
                fastpm_packing_plan_pack(plan, p, j + p->np, (char*) send_buffer + ipar * elsize);
            }
        }

        size_t Nsendsum;
        size_t Nsendallsum;
        MPI_Allreduce(&Nsend, &Nsendsum, 1, MPI_LONG, MPI_SUM, comm);
        MPI_Allreduce(&Nsend_all, &Nsendallsum, 1, MPI_LONG, MPI_SUM, comm);
        fastpm_info("Decomposition iter %d,  exchange of %td particles; need %td", iter, Nsendsum, Nsendallsum);

        MPI_Datatype PTYPE;
        MPI_Type_contiguous(elsize, MPI_BYTE, &PTYPE);
        MPI_Type_commit(&PTYPE);

        MPI_Alltoallv_sparse(
                send_buffer, sendcount, sendoffset, PTYPE,
                recv_buffer, recvcount, recvoffset, PTYPE,
                comm);

        MPI_Type_free(&PTYPE);

        ipar = 0;
        for(i = 0; i < NTask; i ++) {
            for(j = recvoffset[i]; j < recvoffset[i] + recvcount[i]; j++, ipar ++) {
                fastpm_packing_plan_unpack(plan, p, j + p->np, (char*) recv_buffer + ipar * elsize);
            }
        }


        fastpm_memory_free(p->mem, recv_buffer);
        fastpm_memory_free(p->mem, send_buffer);

        p->np += Nrecv;
        iter++;
        if(MPIU_Any(comm, incomplete)) continue;
        else break;
    }
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

int fastpm_store_has_q(FastPMStore * p)
{
    return p->meta._q_size != 0;
}

void
fastpm_store_get_q_from_id(FastPMStore * p, uint64_t id, double q[3])
{
    ptrdiff_t pabs[3];

    int d;
    id = id % p->meta._q_size;
    for(d = 0; d < 3; d++) {
        pabs[d] = id / p->meta._q_strides[d];
        id -= pabs[d] * p->meta._q_strides[d];
    }

    for(d = 0; d < 3; d++) {
        q[d] = pabs[d] * p->meta._q_scale[d];
        q[d] += p->meta._q_shift[d];
    }
}//if all 0 drop

void
fastpm_store_get_iq_from_id(FastPMStore * p, uint64_t id, ptrdiff_t pabs[3])
{
    /* integer version of the initial position q
    i.e. the lattice coordinates of the cells. */
    int d;
    for(d = 0; d < 3; d++) {
        pabs[d] = id / p->meta._q_strides[d];
        id -= pabs[d] * p->meta._q_strides[d];
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
            p->meta._q_shift[d] = shift[d];
        else
            p->meta._q_shift[d] = 0;

        p->meta._q_scale[d] = pm->BoxSize[d] / Nc[d];
    }

    p->meta._q_size = Nc[0] * Nc[1] * Nc[2];
    p->meta._q_strides[0] = Nc[1] * Nc[2];
    p->meta._q_strides[1] = Nc[2];
    p->meta._q_strides[2] = 1;

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

            uint64_t id = pabs[2] * p->meta._q_strides[2] +
                          pabs[1] * p->meta._q_strides[1] +
                          pabs[0] * p->meta._q_strides[0] ;

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
    p->meta.a_x = p->meta.a_v = 0.;
}

void
fastpm_store_summary(FastPMStore * p,
        FastPMColumnTags attribute,
        MPI_Comm comm,
        const char * fmt,
        ...)
{
    va_list va;
    va_start(va, fmt);

    int ci = fastpm_store_find_column_id(p, attribute);
    size_t nmemb = p->_column_info[ci].nmemb;

    double rmin[nmemb], rmax[nmemb], rsum1[nmemb], rsum2[nmemb];

    if (NULL == p->_column_info[ci].to_double) {
        fastpm_raise(-1, "Column %s didnot set to_double virtual function\n",
            p->_column_info[ci].name);
    }
    int d;
    for(d = 0; d < nmemb; d ++) {
        rmin[d] = 1e20;
        rmax[d] = -1e20;
        rsum1[d] = 0;
        rsum2[d] = 0;
    }
    ptrdiff_t i;

#pragma omp parallel
    {
        double tmin[nmemb], tmax[nmemb], tsum1[nmemb], tsum2[nmemb];
        int d;
        for(d = 0; d < nmemb; d ++) {
            tmin[d] = 1e20;
            tmax[d] = -1e20;
            tsum1[d] = 0;
            tsum2[d] = 0;
        }
        #pragma omp for
        for(i = 0; i < p->np; i ++) {
            int d;
            for(d =0; d < nmemb; d++) {
                double value = p->_column_info[ci].to_double(p, i, ci, d);
                tsum1[d] += value;
                tsum2[d] += value * value;
                tmin[d] = fmin(tmin[d], value);
                tmax[d] = fmax(tmax[d], value);
            } 
        }
        #pragma omp critical
        {
            int d;
            for(d =0; d < nmemb; d++) {
                rsum1[d] += tsum1[d];
                rsum2[d] += tsum2[d];
                rmin[d] = fmin(rmin[d], tmin[d]);
                rmax[d] = fmax(rmax[d], tmax[d]);
            }
        }
    }
    uint64_t Ntot = p->np;

    MPI_Allreduce(MPI_IN_PLACE, rsum1, nmemb, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, rsum2, nmemb, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, rmin, nmemb, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, rmax, nmemb, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ntot,   1, MPI_LONG,  MPI_SUM, comm);

    for(i = 0; i < strlen(fmt); i ++) {
        void * r = va_arg(va, void *);
        double * dr = (double * ) r;

        for(d = 0; d < 3; d++) {
            switch(fmt[i]) {
                case '-':
                    dr[d] = rsum1[d] / Ntot;
                    break;
                case '<':
                    dr[d] = rmin[d];
                    break;
                case '>':
                    dr[d] = rmax[d];
                    break;
                case 's':
                    dr[d] = sqrt(rsum2[d] / Ntot - pow(rsum1[d] / Ntot, 2));
                    break;
                case 'S':
                    dr[d] = sqrt(1.0 * Ntot / (Ntot - 1.)) * sqrt(rsum2[d] / Ntot - pow(rsum1[d] / Ntot, 2));
                    break;
                case 'v':
                    dr[d] = (rsum2[d] / Ntot - pow(rsum1[d] / Ntot, 2));
                    break;
                case 'V':
                    dr[d] = (1.0 * Ntot / (Ntot - 1.)) * (rsum2[d] / Ntot - pow(rsum1[d] / Ntot, 2));
                    break;
                default:
                    fastpm_raise(-1, "Unknown format str. Use '<->sSvV'\n");
            }
        }
    }
    va_end(va);
}

void
fastpm_store_steal(FastPMStore * p, FastPMStore * po, FastPMColumnTags attributes)
{
    int c;
    for(c = 0; c < 32; c ++) {
        if (!(p->_column_info[c].attribute & attributes)) continue;
        po->columns[c] = p->columns[c];
    }

    po->np = p->np;
    po->meta = p->meta;
}

static void
_fastpm_store_copy(FastPMStore * p, ptrdiff_t start, FastPMStore * po, ptrdiff_t offset, size_t ncopy)
{
    if(ncopy + start > p->np) {
        fastpm_raise(-1, "Copy out of bounds from source FastPMStore: asking for %td but has %td\n", ncopy + start, p->np);
    }
    if(ncopy + offset > po->np_upper) {
        fastpm_raise(-1, "Not enough storage in target FastPMStore: asking for %td but has %td\n", ncopy + offset, po->np_upper);
    }

    int c;
    for(c = 0; c < 32; c ++) {
        if(!po->columns[c]) continue;
        size_t elsize = po->_column_info[c].elsize;

        memcpy(po->columns[c] + offset * elsize,
                p->columns[c] + start * elsize, elsize * ncopy);

    }
    po->np = offset + ncopy;
    po->meta = p->meta;
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

/* extends p by extra. */
void
fastpm_store_extend(FastPMStore * p, FastPMStore * extra)
{
    _fastpm_store_copy(extra, 0, p, p->np, extra->np);
}

void
fastpm_store_fill_subsample_mask(FastPMStore * p,
        double fraction,
        FastPMParticleMaskType * mask,
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

    memset(mask, 0, p->np * sizeof(mask[0]));

    ptrdiff_t i;
    for(i=0; i < p->np; i++) {
        double rand_i = gsl_rng_uniform(random_generator);
        int flag = fraction > 1 || rand_i <= fraction;
        mask[i] = flag;
    }

    gsl_rng_free(random_generator);
}

void
fastpm_store_fill_subsample_mask_every_dim(FastPMStore * p,
                                              int every, /* take 1 every 'every' per dimension */
                                              FastPMParticleMaskType * mask)
{
    if(!fastpm_store_has_q(p)) {
        /* This can be relaxed by using a less strict subsample algorithm, e.g. subsample by a hash of ID. */
        fastpm_raise(-1, "Subsample is not supported if the store does not have meta.q.");
    }
    /* UNUSED */
    memset(mask, 0, p->np * sizeof(mask[0]));

    ptrdiff_t i, d;
    for(i = 0; i < p->np; i++) {
        ptrdiff_t pabs[3];   //pabs is the iq index. move out of loop?
        uint64_t id = p->id[i];
        fastpm_store_get_iq_from_id(p, id, pabs);

        int flag = 1;
        for(d = 0; d < 3; d++) {
            flag *= !(pabs[d] % every);
        }
        mask[i] = flag;
    }
}

/*
 * Create a subsample, keeping only those with mask == True;
 *
 * if po is NULL, only return number of items.
 * */
size_t
fastpm_store_subsample(FastPMStore * p, FastPMParticleMaskType * mask, FastPMStore * po)
{
    ptrdiff_t i;
    ptrdiff_t j;

    j = 0;
    for(i = 0; i < p->np; i ++) {
        if(!mask[i]) continue;
        /* just counting */
        if(po == NULL) { j++; continue; }
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

    if(po) {
        po->np = j;
        po->meta = p->meta;   ///????
    }
    return j;
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
