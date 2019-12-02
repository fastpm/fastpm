
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <mpi.h>

#include "mpsort.h"
#include "internal.h"

#include "mp-mpiu.h"
#include "internal-parallel.h"

static int _mpsort_mpi_options = 0;

/* mpi version of radix sort;
 *
 * each caller provides the distributed array and number of items.
 * the sorted array is returned to the original array pointed to by
 * mybase. (AKA no rebalancing is done.)
 *
 * NOTE: may need an api to return a balanced array!
 *
 * uses the same amount of temporary storage space for communication
 * and local sort. (this will be allocated via malloc)
 *
 *
 * */

static MPI_Datatype MPI_TYPE_PTRDIFF = 0;

struct crmpistruct {
    MPI_Datatype MPI_TYPE_RADIX;
    MPI_Datatype MPI_TYPE_DATA;
    MPI_Comm comm;
    void * mybase;
    void * myoutbase;
    size_t mynmemb;
    size_t nmemb;
    size_t myoutnmemb;
    size_t outnmemb;
    int NTask;
    int ThisTask;
};

static void
_setup_mpsort_mpi(struct crmpistruct * o,
                  struct crstruct * d,
                  void * myoutbase, size_t myoutnmemb,
                  MPI_Comm comm, const int line, const char * file)
{

    o->comm = comm;

    MPI_Comm_size(comm, &o->NTask);
    MPI_Comm_rank(comm, &o->ThisTask);

    o->mybase = d->base;
    o->mynmemb = d->nmemb;
    o->myoutbase = myoutbase;
    o->myoutnmemb = myoutnmemb;

    MPI_Allreduce(&o->mynmemb, &o->nmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, comm);
    MPI_Allreduce(&o->myoutnmemb, &o->outnmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, comm);

    if(o->outnmemb != o->nmemb) {
        if(o->ThisTask == 0) {
            fprintf(stderr, "MPSort: total number of items in the item does not match the input %ld != %ld. "
                            "Caller site: %s:%d\n",
                            o->outnmemb, o->nmemb, file, line);
            MPI_Abort(comm, -1);
        }
    }

    MPI_Type_contiguous(d->rsize, MPI_BYTE, &o->MPI_TYPE_RADIX);
    MPI_Type_commit(&o->MPI_TYPE_RADIX);

    MPI_Type_contiguous(d->size, MPI_BYTE, &o->MPI_TYPE_DATA);
    MPI_Type_commit(&o->MPI_TYPE_DATA);

}
static void _destroy_mpsort_mpi(struct crmpistruct * o) {
    MPI_Type_free(&o->MPI_TYPE_RADIX);
    MPI_Type_free(&o->MPI_TYPE_DATA);
}

static void _find_Pmax_Pmin_C(void * mybase, size_t mynmemb, size_t nmemb,
        size_t myoutnmemb,
        unsigned char * Pmax, unsigned char * Pmin,
        ptrdiff_t * C,
        struct crstruct * d,
        struct crmpistruct * o);

static int _solve_for_layout_mpi (
        int NTask,
        ptrdiff_t * C,
        ptrdiff_t * myT_CLT,
        ptrdiff_t * myT_CLE,
        ptrdiff_t * myT_C,
        MPI_Comm comm);

static struct TIMER {
    double time;
    char name[20];
} _TIMERS[512];

void mpsort_mpi_report_last_run() {
    struct TIMER * tmr = _TIMERS;
    double last = tmr->time;
    tmr ++;
    while(0 != strcmp(tmr->name, "END")) {
        printf("%s: %g\n", tmr->name, tmr->time - last);
        last =tmr->time;
        tmr ++;
    }
}
int mpsort_mpi_find_ntimers(struct TIMER * tmr) {
    int n = 0;
    while(0 != strcmp(tmr->name, "END")) {
        tmr ++;
        n++;
    }
    return n;
}

void
mpsort_mpi_impl (void * mybase, size_t mynmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg,
        MPI_Comm comm,
        const int line,
        const char * file)
{

    mpsort_mpi_newarray_impl(mybase, mynmemb,
        mybase, mynmemb,
        size, radix, rsize, arg, comm, line, file);
}

static int
mpsort_mpi_histogram_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr,
    const int line, const char * file);

static uint64_t
checksum(void * base, ptrdiff_t nbytes, MPI_Comm comm)
{
    uint64_t sum = 0;
    char * ptr = base;
    ptrdiff_t i = 0;
    for(i = 0; i < nbytes; i ++) {
        sum += ptr[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_LONG, MPI_SUM, comm);
    return sum;
}

void
mpsort_mpi_newarray_impl (void * mybase, size_t mynmemb,
        void * myoutbase, size_t myoutnmemb,
        size_t elsize,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg,
        MPI_Comm comm,
        const int line,
        const char * file)
{

    if(MPI_TYPE_PTRDIFF == 0) {
        if(sizeof(ptrdiff_t) == sizeof(int)) {
            MPI_TYPE_PTRDIFF = MPI_INT;
        }
        else if(sizeof(ptrdiff_t) == sizeof(long)) {
            MPI_TYPE_PTRDIFF = MPI_LONG;
        }
        else if(sizeof(ptrdiff_t) == sizeof(long long)) {
            MPI_TYPE_PTRDIFF = MPI_LONG_LONG;
        }
        else {
            fprintf(stderr, "MPSort: sizeof(ptrdiff) = %lu, not recognised\n", sizeof(ptrdiff_t));
            MPI_Abort(comm, -1);
        }
    }

    struct TIMER * tmr = _TIMERS;

    MPIU_Segmenter segmenter[1];

    uint64_t sum1 = checksum(mybase, elsize * mynmemb, comm);

    int NTask;
    int ThisTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    if(elsize > 8 && elsize % 8 != 0) {
        if(ThisTask == 0) {
            fprintf(stderr, "MPSort: element size is large (%ld) but not aligned to 8 bytes. "
                            "This is known to frequently trigger MPI bugs. "
                            "Caller site: %s:%d\n",
                            elsize, file, line);
        }
    }
    if(rsize > 8 && rsize % 8 != 0) {
        if(ThisTask == 0) {
            fprintf(stderr, "MPSort: radix size is large (%ld) but not aligned to 8 bytes. "
                            "This is known to frequently trigger MPI bugs. "
                            "Caller site: %s:%d\n",
                            rsize, file, line);
        }
    }

    size_t sizes[NTask];
    size_t outsizes[NTask];
    size_t myoffset;
    size_t myoutoffset;
    /* comput sizes for both input and output to avoid any of them getting too imbalanced */
    size_t totalsize = MPIU_Segmenter_collect_sizes(mynmemb, sizes, &myoffset, comm);
    size_t totalsizeout = MPIU_Segmenter_collect_sizes(myoutnmemb, outsizes, &myoutoffset, comm);

    if(totalsize != totalsizeout) {
        if(ThisTask == 0) {
            fprintf(stderr, "Input and output size mismatch: %td (in) != %td (out)"
                            "Caller site: %s:%d\n",
                            totalsize, totalsizeout, file, line);
        }
        MPI_Abort(comm, -1);
    }

    size_t avgsegsize = NTask; /* combine very small ranks to segments */
    if (avgsegsize * elsize > 4 * 1024 * 1024) {
        /* do not use more than 4MB in a segment */
        avgsegsize = 4 * 1024 * 1024 / elsize;
    }
    if(mpsort_mpi_has_options(MPSORT_REQUIRE_GATHER_SORT)) {
        if(ThisTask == 0) {
            fprintf(stderr, "MPSort: gathering all data to a single rank for sorting due to MPSORT_REQUIRE_GATHER_SORT. "
                            "Total number of items is %ld. "
                            "Caller site: %s:%d\n",
                            totalsize, file, line);
        }
        avgsegsize = totalsize;
    }

    if(mpsort_mpi_has_options(MPSORT_DISABLE_GATHER_SORT)) {
        avgsegsize = 0;
        if(ThisTask == 0) {
            fprintf(stderr, "MPSort: disable gathering data into larger chunks due to MPSORT_DISABLE_GATHER_SORT. "
                            "Caller site: %s:%d\n",
                            file, line);
        }
    }

    /* use as many groups as possible (some will be empty) but at most 1 segment per group */
    MPIU_Segmenter_init(segmenter, sizes, outsizes, avgsegsize, NTask, comm);

    /* group comm == seg comm */

    void * mysegmentbase = NULL;
    void * myoutsegmentbase = NULL;
    size_t mysegmentnmemb;
    size_t myoutsegmentnmemb;

    int groupsize;
    int grouprank;
    MPI_Comm_size(segmenter->Group, &groupsize);
    MPI_Comm_rank(segmenter->Group, &grouprank);

    MPI_Allreduce(&mynmemb, &mysegmentnmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, segmenter->Group);
    MPI_Allreduce(&myoutnmemb, &myoutsegmentnmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, segmenter->Group);

    if (groupsize > 1) {
        if(grouprank == segmenter->group_leader_rank) {
            mysegmentbase = MPIU_Malloc("mysegment", elsize, mysegmentnmemb);
            myoutsegmentbase = MPIU_Malloc("outsegment", elsize, myoutsegmentnmemb);
        }
        MPIU_Gather(segmenter->Group, segmenter->group_leader_rank, mybase, mysegmentbase, mynmemb, elsize, NULL);
    } else {
        mysegmentbase = mybase;
        myoutsegmentbase = myoutbase;
    }

    /* only do sorting on the group leaders for each segment */
    if(segmenter->is_group_leader) {

        struct crstruct d;
        struct crmpistruct o;

        _setup_radix_sort(&d, mysegmentbase, mysegmentnmemb, elsize, radix, rsize, arg);

        _setup_mpsort_mpi(&o, &d, myoutsegmentbase, myoutsegmentnmemb, segmenter->Leaders, line, file);

        mpsort_mpi_histogram_sort(d, o, tmr, line, file);

        _destroy_mpsort_mpi(&o);
    }

    if(groupsize > 1) {
        MPIU_Scatter(segmenter->Group, segmenter->group_leader_rank, myoutsegmentbase, myoutbase, myoutnmemb, elsize, NULL);
    }

    {
        int ntmr;
        if(segmenter->is_group_leader)
            ntmr = (mpsort_mpi_find_ntimers(tmr) + 1);

        MPI_Bcast(&ntmr, 1, MPI_INT, segmenter->group_leader_rank, segmenter->Group);
        MPI_Bcast(tmr, sizeof(tmr[0]) * ntmr, MPI_BYTE, segmenter->group_leader_rank, segmenter->Group);
    }

    if(grouprank == segmenter->group_leader_rank) {
        if(mysegmentbase != mybase)
            MPIU_Free(mysegmentbase);
        if(myoutsegmentbase != myoutbase)
            MPIU_Free(myoutsegmentbase);
    }

    MPIU_Segmenter_destroy(segmenter);

    uint64_t sum2 = checksum(myoutbase, elsize * myoutnmemb, comm);
    if (sum1 != sum2) {
        fprintf(stderr, "MPSort: Data changed after sorting; checksum mismatch. "
                        "Caller site: %s:%d\n",
                        file, line);
        MPI_Abort(comm, -1);
    }
}

int
mpsort_mpi_histogram_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr,
        const int line, const char * file)
{

    unsigned char Pmax[d.rsize];
    unsigned char Pmin[d.rsize];

    unsigned char P[d.rsize * (o.NTask - 1)];

    ptrdiff_t C[o.NTask + 1];  /* desired counts */

    ptrdiff_t myCLT[o.NTask + 1]; /* counts of less than P */
    ptrdiff_t CLT[o.NTask + 1];

    ptrdiff_t myCLE[o.NTask + 1]; /* counts of less than or equal to P */
    ptrdiff_t CLE[o.NTask + 1];

    int SendCount[o.NTask];
    int SendDispl[o.NTask];
    int RecvCount[o.NTask];
    int RecvDispl[o.NTask];

    ptrdiff_t myT_CLT[o.NTask];
    ptrdiff_t myT_CLE[o.NTask];
    ptrdiff_t myT_C[o.NTask];
    ptrdiff_t myC[o.NTask + 1];

    int iter = 0;
    int done = 0;
    char * buffer;
    int i;

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "START"), tmr++);

    /* and sort the local array */
    radix_sort(d.base, d.nmemb, d.size, d.radix, d.rsize, d.arg);

    MPI_Barrier(o.comm);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "FirstSort"), tmr++);

    _find_Pmax_Pmin_C(o.mybase, o.mynmemb, o.nmemb, o.myoutnmemb, Pmax, Pmin, C, &d, &o);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "PmaxPmin"), tmr++);

    memset(P, 0, d.rsize * (o.NTask -1));

    struct piter pi;

    piter_init(&pi, Pmin, Pmax, o.NTask - 1, &d);

    while(!done) {
        iter ++;
        piter_bisect(&pi, P);

        _histogram(P, o.NTask - 1, o.mybase, o.mynmemb, myCLT, myCLE, &d);

        MPI_Allreduce(myCLT, CLT, o.NTask + 1,
                MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);
        MPI_Allreduce(myCLE, CLE, o.NTask + 1,
                MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);

        (iter>10?tmr--:0, tmr->time = MPI_Wtime(), sprintf(tmr->name, "bisect%04d", iter), tmr++);

        piter_accept(&pi, P, C, CLT, CLE);
#if 0
        {
            int k;
            for(k = 0; k < o.NTask; k ++) {
                MPI_Barrier(o.comm);
                int i;
                if(o.ThisTask != k) continue;

                printf("P (%d): PMin %d PMax %d P ",
                        o.ThisTask,
                        *(int*) Pmin,
                        *(int*) Pmax
                        );
                for(i = 0; i < o.NTask - 1; i ++) {
                    printf(" %d ", ((int*) P) [i]);
                }
                printf("\n");

                printf("C (%d): ", o.ThisTask);
                for(i = 0; i < o.NTask + 1; i ++) {
                    printf("%d ", C[i]);
                }
                printf("\n");
                printf("CLT (%d): ", o.ThisTask);
                for(i = 0; i < o.NTask + 1; i ++) {
                    printf("%d ", CLT[i]);
                }
                printf("\n");
                printf("CLE (%d): ", o.ThisTask);
                for(i = 0; i < o.NTask + 1; i ++) {
                    printf("%d ", CLE[i]);
                }
                printf("\n");

            }
        }
#endif
        done = piter_all_done(&pi);
    }

    piter_destroy(&pi);

    _histogram(P, o.NTask - 1, o.mybase, o.mynmemb, myCLT, myCLE, &d);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "findP"), tmr++);

    /* transpose the matrix, could have been done with a new datatype */
    /*
    MPI_Alltoall(myCLT, 1, MPI_TYPE_PTRDIFF,
            myT_CLT, 1, MPI_TYPE_PTRDIFF, o.comm);
    */
    MPI_Alltoall(myCLT + 1, 1, MPI_TYPE_PTRDIFF,
            myT_CLT, 1, MPI_TYPE_PTRDIFF, o.comm);

    /*MPI_Alltoall(myCLE, 1, MPI_TYPE_PTRDIFF,
            myT_CLE, 1, MPI_TYPE_PTRDIFF, o.comm); */
    MPI_Alltoall(myCLE + 1, 1, MPI_TYPE_PTRDIFF,
            myT_CLE, 1, MPI_TYPE_PTRDIFF, o.comm);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "LayDistr"), tmr++);

    _solve_for_layout_mpi(o.NTask, C, myT_CLT, myT_CLE, myT_C, o.comm);

    myC[0] = 0;
    MPI_Alltoall(myT_C, 1, MPI_TYPE_PTRDIFF,
            myC + 1, 1, MPI_TYPE_PTRDIFF, o.comm);
#if 0
    for(i = 0;i < o.NTask; i ++) {
        int j;
        MPI_Barrier(o.comm);
        if(o.ThisTask != i) continue;
        for(j = 0; j < o.NTask + 1; j ++) {
            printf("%d %d %d, ",
                    myCLT[j],
                    myC[j],
                    myCLE[j]);
        }
        printf("\n");

    }
#endif

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "LaySolve"), tmr++);

    for(i = 0; i < o.NTask; i ++) {
        SendCount[i] = myC[i + 1] - myC[i];
    }

    MPI_Alltoall(SendCount, 1, MPI_INT,
            RecvCount, 1, MPI_INT, o.comm);

    SendDispl[0] = 0;
    RecvDispl[0] = 0;
    size_t totrecv = RecvCount[0];
    for(i = 1; i < o.NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
        if(SendDispl[i] != myC[i]) {
            fprintf(stderr, "SendDispl: error. "
                        "Caller site: %s:%d\n",
                        file, line);
            MPI_Abort(o.comm, -1);
        }
        totrecv += RecvCount[i];
    }
    if(totrecv != o.myoutnmemb) {
        fprintf(stderr, "totrecv = %td, mismatch with %td. "
                        "Caller site: %s:%d\n",
                        totrecv, o.myoutnmemb,
                        file, line);
        MPI_Abort(o.comm, -1);
    }
#if 0
    {
        int k;
        for(k = 0; k < o.NTask; k ++) {
            MPI_Barrier(o.comm);

            if(o.ThisTask != k) continue;

            printf("P (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask - 1; i ++) {
                printf("%d ", ((int*) P) [i]);
            }
            printf("\n");

            printf("C (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", C[i]);
            }
            printf("\n");
            printf("CLT (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", CLT[i]);
            }
            printf("\n");
            printf("CLE (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", CLE[i]);
            }
            printf("\n");

            printf("MyC (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", myC[i]);
            }
            printf("\n");
            printf("MyCLT (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", myCLT[i]);
            }
            printf("\n");

            printf("MyCLE (%d): ", o.ThisTask);
            for(i = 0; i < o.NTask + 1; i ++) {
                printf("%d ", myCLE[i]);
            }
            printf("\n");

            printf("Send Count(%d): ", o.ThisTask);
            for(i = 0; i < o.NTask; i ++) {
                printf("%d ", SendCount[i]);
            }
            printf("\n");
            printf("My data(%d): ", o.ThisTask);
            for(i = 0; i < mynmemb; i ++) {
                printf("%d ", ((int*) mybase)[i]);
            }
            printf("\n");
        }
    }
#endif
    if(o.myoutbase == o.mybase)
        buffer = MPIU_Malloc("buffer", d.size, o.myoutnmemb);
    else
        buffer = o.myoutbase;

    enum MPIU_AlltoallvSparsePolicy policy = AUTO;
    if (mpsort_mpi_has_options(MPSORT_DISABLE_SPARSE_ALLTOALLV)) {
        policy = DISABLED;
    }
    if (mpsort_mpi_has_options(MPSORT_REQUIRE_SPARSE_ALLTOALLV)) {
        policy = REQUIRED;
    }

    MPIU_Alltoallv(
            o.mybase, SendCount, SendDispl, o.MPI_TYPE_DATA,
            buffer, RecvCount, RecvDispl, o.MPI_TYPE_DATA,
            o.comm, policy);

    if(o.myoutbase == o.mybase) {
        memcpy(o.myoutbase, buffer, o.myoutnmemb * d.size);
        MPIU_Free(buffer);
    }

    MPI_Barrier(o.comm);
    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "Exchange"), tmr++);

    radix_sort(o.myoutbase, o.myoutnmemb, d.size, d.radix, d.rsize, d.arg);

    MPI_Barrier(o.comm);
    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "SecondSort"), tmr++);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "END"), tmr++);
    return 0;
}

static void _find_Pmax_Pmin_C(void * mybase, size_t mynmemb,
        size_t nmemb,
        size_t myoutnmemb,
        unsigned char * Pmax, unsigned char * Pmin,
        ptrdiff_t * C,
        struct crstruct * d,
        struct crmpistruct * o) {
    memset(Pmax, 0, d->rsize);
    memset(Pmin, -1, d->rsize);

    unsigned char myPmax[d->rsize];
    unsigned char myPmin[d->rsize];

    size_t eachnmemb[o->NTask];
    size_t eachoutnmemb[o->NTask];
    unsigned char eachPmax[d->rsize * o->NTask];
    unsigned char eachPmin[d->rsize * o->NTask];
    int i;

    if(mynmemb > 0) {
        d->radix((char*) mybase + (mynmemb - 1) * d->size, myPmax, d->arg);
        d->radix(mybase, myPmin, d->arg);
    } else {
        memset(myPmin, 0, d->rsize);
        memset(myPmax, 0, d->rsize);
    }

    MPI_Allgather(&mynmemb, 1, MPI_TYPE_PTRDIFF,
            eachnmemb, 1, MPI_TYPE_PTRDIFF, o->comm);
    MPI_Allgather(&myoutnmemb, 1, MPI_TYPE_PTRDIFF,
            eachoutnmemb, 1, MPI_TYPE_PTRDIFF, o->comm);
    MPI_Allgather(myPmax, 1, o->MPI_TYPE_RADIX,
            eachPmax, 1, o->MPI_TYPE_RADIX, o->comm);
    MPI_Allgather(myPmin, 1, o->MPI_TYPE_RADIX,
            eachPmin, 1, o->MPI_TYPE_RADIX, o->comm);


    C[0] = 0;
    for(i = 0; i < o->NTask; i ++) {
        C[i + 1] = C[i] + eachoutnmemb[i];
        /* skip the rank, since it has no data to contribute to the reduction */
        if(eachnmemb[i] == 0) continue;

        if(d->compar(eachPmax + i * d->rsize, Pmax, d->rsize) > 0) {
            memcpy(Pmax, eachPmax + i * d->rsize, d->rsize);
        }
        if(d->compar(eachPmin + i * d->rsize, Pmin, d->rsize) < 0) {
            memcpy(Pmin, eachPmin + i * d->rsize, d->rsize);
        }
    }
    /* no rank contributed to the reduction, set min, max to a sane value, 0 */
    if(nmemb == 0) {
        memset(Pmin, 0, d->rsize);
        memset(Pmax, 0, d->rsize);
    }
}

static int
_solve_for_layout_mpi (
        int NTask,
        ptrdiff_t * C,
        ptrdiff_t * myT_CLT,
        ptrdiff_t * myT_CLE,
        ptrdiff_t * myT_C,
        MPI_Comm comm) {
    int i, j;
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    /* first assume we just send according to myT_CLT */
    for(i = 0; i < NTask; i ++) {
        myT_C[i] = myT_CLT[i];
    }

    /* Solve for each receiving task i
     *
     * this solves for GL_C[..., i + 1], which depends on GL_C[..., i]
     *
     * and we have GL_C[..., 0] == 0 by definition.
     *
     * this cannot be done in parallel wrt i because of the dependency.
     *
     *  a solution is guaranteed because GL_CLE and GL_CLT
     *  brackes the total counts C (we've found it with the
     *  iterative counting.
     *
     * */

    ptrdiff_t sure = 0;

    /* how many will I surely receive? */
    for(j = 0; j < NTask; j ++) {
        ptrdiff_t recvcount = myT_C[j];
        sure += recvcount;
    }
    /* let's see if we have enough */
    ptrdiff_t deficit = C[ThisTask + 1] - sure;

    for(j = 0; j < NTask; j ++) {
        /* deficit solved */
        if(deficit == 0) break;
        if(deficit < 0) {
            fprintf(stderr, "serious bug: more items than there should be: deficit=%ld\n", deficit);
            abort();
        }
        /* how much task j can supply ? */
        ptrdiff_t supply = myT_CLE[j] - myT_C[j];
        if(supply < 0) {
            fprintf(stderr, "serious bug: less items than there should be: supply =%ld\n", supply);
            abort();
        }
        if(supply <= deficit) {
            myT_C[j] += supply;
            deficit -= supply;
        } else {
            myT_C[j] += deficit;
            deficit = 0;
        }
    }

    return 0;
}



static void _mpsort_mpi_parse_env()
{
    static int _mpsort_env_parsed = 0;
    if(_mpsort_env_parsed) return;

    _mpsort_env_parsed = 1;
    if(getenv("MPSORT_DISABLE_SPARSE_ALLTOALLV"))
        mpsort_mpi_set_options(MPSORT_DISABLE_SPARSE_ALLTOALLV);
    if(getenv("MPSORT_DISABLE_GATHER_SORT"))
        mpsort_mpi_set_options(MPSORT_DISABLE_GATHER_SORT);
    if(getenv("MPSORT_REQUIRE_GATHER_SORT "))
        mpsort_mpi_set_options(MPSORT_REQUIRE_GATHER_SORT );
    if(getenv("MPSORT_REQUIRE_SPARSE_ALLTOALLV"))
        mpsort_mpi_set_options(MPSORT_REQUIRE_SPARSE_ALLTOALLV);
}

void
mpsort_mpi_set_options(int options)
{
    _mpsort_mpi_parse_env();
    _mpsort_mpi_options |= options;
}

int
mpsort_mpi_has_options(int options)
{
    _mpsort_mpi_parse_env();
    return _mpsort_mpi_options & options;
}

void
mpsort_mpi_unset_options(int options)
{
    _mpsort_mpi_parse_env();
    _mpsort_mpi_options &= ~options;

}
