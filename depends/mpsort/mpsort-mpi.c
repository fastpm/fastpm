
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <mpi.h>

#include "mpsort.h"

#include "internal.h"

#include "internal-parallel.h"

static int _mpsort_mpi_options = 0;

static int MPI_Alltoallv_smart(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);

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
                  MPI_Comm comm)
{

    o->comm = comm;

    MPI_Comm_size(comm, &o->NTask);
    MPI_Comm_rank(comm, &o->ThisTask);

    if(MPI_TYPE_PTRDIFF == 0) {
        if(sizeof(ptrdiff_t) == sizeof(long long)) {
            MPI_TYPE_PTRDIFF = MPI_LONG_LONG;
        }
        if(sizeof(ptrdiff_t) == sizeof(long)) {
            MPI_TYPE_PTRDIFF = MPI_LONG;
        }
        if(sizeof(ptrdiff_t) == sizeof(int)) {
            MPI_TYPE_PTRDIFF = MPI_INT;
        }
    }

    o->mybase = d->base;
    o->mynmemb = d->nmemb;
    o->myoutbase = myoutbase;
    o->myoutnmemb = myoutnmemb;

    MPI_Allreduce(&o->mynmemb, &o->nmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, comm);
    MPI_Allreduce(&o->myoutnmemb, &o->outnmemb, 1, MPI_TYPE_PTRDIFF, MPI_SUM, comm);

    if(o->outnmemb != o->nmemb) {
        fprintf(stderr, "total number of items in the item does not match the input %ld != %ld\n",
                o->outnmemb, o->nmemb);
        abort();
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

static void _find_Pmax_Pmin_C(void * mybase, size_t mynmemb, 
        size_t myoutnmemb,
        char * Pmax, char * Pmin, 
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
} _TIMERS[20];

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

void
mpsort_mpi (void * mybase, size_t mynmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg, 
        MPI_Comm comm)
{

    mpsort_mpi_newarray(mybase, mynmemb,
        mybase, mynmemb, 
        size, radix, rsize, arg, comm);
}

int
mpsort_mpi_gather_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr);

int
mpsort_mpi_histogram_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr);

void
mpsort_mpi_newarray (void * mybase, size_t mynmemb, 
        void * myoutbase, size_t myoutnmemb,
        size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg,
        MPI_Comm comm)
{

    struct TIMER * tmr = _TIMERS;

    struct crstruct d;
    struct crmpistruct o;

    _setup_radix_sort(&d, mybase, mynmemb, size, radix, rsize, arg);
    _setup_mpsort_mpi(&o, &d, myoutbase, myoutnmemb, comm);

    if(o.nmemb == 0) goto exec_empty_array;

    if( mpsort_mpi_has_options(MPSORT_REQUIRE_GATHER_SORT) ||
      (!mpsort_mpi_has_options(MPSORT_DISABLE_GATHER_SORT)
        && o.nmemb <= o.NTask && d.size <= 32)) {
        mpsort_mpi_gather_sort(d, o, tmr);
    } else {
        mpsort_mpi_histogram_sort(d, o, tmr);
    }
exec_empty_array:
    _destroy_mpsort_mpi(&o);

}


int
mpsort_mpi_gather_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr)
{
    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "START"), tmr++);

    int recvcounts[o.NTask + 1];
    int rdispls[o.NTask + 1];
    int sendcounts[o.NTask + 1];
    int sdispls[o.NTask + 1];

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "START"), tmr++);

    MPI_Gather(&o.mynmemb, 1, MPI_INT,
                 recvcounts, 1, MPI_INT, 0, o.comm);

    MPI_Gather(&o.myoutnmemb, 1, MPI_INT,
                 sendcounts, 1, MPI_INT, 0, o.comm);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "layout"), tmr++);

    rdispls[0] = 0;
    sdispls[0] = 0;
    int i;
    for (i = 1; i < o.NTask; i++) {
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
        sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
    }

    void * buffer = malloc(d.size * o.nmemb);

    MPI_Gatherv(o.mybase, o.mynmemb, o.MPI_TYPE_DATA, buffer, recvcounts, rdispls, o.MPI_TYPE_DATA, 0, o.comm);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "Gather"), tmr++);
    /* and sort the local array */
    if (o.ThisTask == 0) {
        radix_sort(buffer, o.nmemb, d.size, d.radix, d.rsize, d.arg);
    }
    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "LocalSort"), tmr++);

    MPI_Scatterv(buffer, sendcounts, sdispls, o.MPI_TYPE_DATA, o.myoutbase, o.myoutnmemb, o.MPI_TYPE_DATA, 0, o.comm);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "Scatter"), tmr++);

    free(buffer);
    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "END"), tmr++);
    return 0;
}

int
mpsort_mpi_histogram_sort(struct crstruct d, struct crmpistruct o, struct TIMER * tmr)
{

    char Pmax[d.rsize];
    char Pmin[d.rsize];

    char P[d.rsize * (o.NTask - 1)];

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

    _find_Pmax_Pmin_C(o.mybase, o.mynmemb, o.myoutnmemb, Pmax, Pmin, C, &d, &o);

    (tmr->time = MPI_Wtime(), strcpy(tmr->name, "PmaxPmin"), tmr++);

    memset(P, 0, d.rsize * (o.NTask -1));

    struct piter pi;

    piter_init(&pi, Pmin, Pmax, o.NTask - 1, &d);

    while(!done) {
        iter ++;
        piter_bisect(&pi, P);

        _histogram(P, o.NTask - 1, o.mybase, o.mynmemb, myCLT, myCLE, &d);

#if MPI_VERSION >= 3
        if (mpsort_mpi_has_options(MPSORT_DISABLE_IALLREDUCE)
        ) {
            MPI_Allreduce(myCLT, CLT, o.NTask + 1, 
                    MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);
            MPI_Allreduce(myCLE, CLE, o.NTask + 1, 
                    MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);
        } else {
            MPI_Request r[2];

            MPI_Iallreduce(myCLT, CLT, o.NTask + 1, 
                    MPI_TYPE_PTRDIFF, MPI_SUM, o.comm, &r[0]);
            MPI_Iallreduce(myCLE, CLE, o.NTask + 1, 
                    MPI_TYPE_PTRDIFF, MPI_SUM, o.comm, &r[1]);

            MPI_Waitall(2, r, MPI_STATUSES_IGNORE);
            /* no free because MPI_Wait frees it automatially. */
        }
#else
        MPI_Allreduce(myCLT, CLT, o.NTask + 1, 
                MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);
        MPI_Allreduce(myCLE, CLE, o.NTask + 1, 
                MPI_TYPE_PTRDIFF, MPI_SUM, o.comm);
#endif

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
            fprintf(stderr, "SendDispl error\n");
            abort();
        }
        totrecv += RecvCount[i];
    }
    if(totrecv != o.myoutnmemb) {
        fprintf(stderr, "totrecv = %td, mismatch with %td\n", totrecv, o.myoutnmemb);
        abort();
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
        buffer = malloc(d.size * o.myoutnmemb);
    else
        buffer = o.myoutbase;

    MPI_Alltoallv_smart(
            o.mybase, SendCount, SendDispl, o.MPI_TYPE_DATA,
            buffer, RecvCount, RecvDispl, o.MPI_TYPE_DATA, 
            o.comm);

    if(o.myoutbase == o.mybase) {
        memcpy(o.myoutbase, buffer, o.myoutnmemb * d.size);
        free(buffer);
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
        size_t myoutnmemb,
        char * Pmax, char * Pmin, 
        ptrdiff_t * C,
        struct crstruct * d,
        struct crmpistruct * o) {
    memset(Pmax, 0, d->rsize);
    memset(Pmin, -1, d->rsize);

    char myPmax[d->rsize];
    char myPmin[d->rsize];

    size_t eachnmemb[o->NTask];
    size_t eachoutnmemb[o->NTask];
    char eachPmax[d->rsize * o->NTask];
    char eachPmin[d->rsize * o->NTask];
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
        if(eachnmemb[i] == 0) continue;

        if(d->compar(eachPmax + i * d->rsize, Pmax, d->rsize) > 0) {
            memcpy(Pmax, eachPmax + i * d->rsize, d->rsize);
        }
        if(d->compar(eachPmin + i * d->rsize, Pmin, d->rsize) < 0) {
            memcpy(Pmin, eachPmin + i * d->rsize, d->rsize);
        }
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



/* The following two functions are taken from MP-Gadget. The hope 
 * is that when the exchange is sparse posting requests is
 * faster than Alltoall on some implementations. */

static int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);

int MPI_Alltoallv_smart(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) 
/* 
 * sdispls, recvcnts rdispls can be NULL,
 *
 * if recvbuf is NULL, returns total number of item required to hold the
 * data.
 * */
{
    int ThisTask;
    int NTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);
    int i;
    int nn = 0;
    int *a_sdispls=NULL, *a_recvcnts=NULL, *a_rdispls=NULL;
    for(i = 0; i < NTask; i ++) {
        if(sendcnts[i] > 0) {
            nn ++;
        }
    }
    if(recvcnts == NULL) {
        a_recvcnts = malloc(sizeof(int) * NTask);
        recvcnts = a_recvcnts;
        MPI_Alltoall(sendcnts, 1, MPI_INT,
                     recvcnts, 1, MPI_INT, comm);
    }
    if(recvbuf == NULL) {
        int totalrecv = 0;
        for(i = 0; i < NTask; i ++) {
            totalrecv += recvcnts[i];
        }
        if(a_recvcnts)
            free(a_recvcnts);
        return totalrecv;
    }
    if(sdispls == NULL) {
        a_sdispls = malloc(sizeof(int) * NTask);
        sdispls = a_sdispls;
        sdispls[0] = 0;
        for (i = 1; i < NTask; i++) {
            sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
        }
    }
    if(rdispls == NULL) {
        a_rdispls = malloc(sizeof(int) * NTask);
        rdispls = a_rdispls;
        rdispls[0] = 0;
        for (i = 1; i < NTask; i++) {
            rdispls[i] = rdispls[i - 1] + recvcnts[i - 1];
        }
    }

    int dense;

    if(mpsort_mpi_has_options(MPSORT_ENABLE_SPARSE_ALLTOALLV)
    ) {
        dense = 1;
    } else {
        dense = nn < NTask * 0.2;
        MPI_Allreduce(MPI_IN_PLACE, &dense, 1, MPI_INT, MPI_SUM, comm);
    }

    if(mpsort_mpi_has_options(MPSORT_REQUIRE_SPARSE_ALLTOALLV)) {
        dense = 0;
    }

    int ret;
    if(dense != 0) {
        ret = MPI_Alltoallv(sendbuf, sendcnts, sdispls,
                    sendtype, recvbuf, 
                    recvcnts, rdispls, recvtype, comm);
    } else {
        ret = MPI_Alltoallv_sparse(sendbuf, sendcnts, sdispls,
                    sendtype, recvbuf, 
                    recvcnts, rdispls, recvtype, comm);
    }

    if(a_rdispls)
        free(a_rdispls);
    if(a_sdispls)
        free(a_sdispls);
    if(a_recvcnts)
        free(a_recvcnts);
    return ret;
}

static int MPI_Alltoallv_sparse(void *sendbuf, int *sendcnts, int *sdispls,
        MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
        int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {

    int ThisTask;
    int NTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);
    int PTask;
    int ngrp;

    for(PTask = 0; NTask > (1 << PTask); PTask++);

    ptrdiff_t lb;
    ptrdiff_t send_elsize;
    ptrdiff_t recv_elsize;

    MPI_Type_get_extent(sendtype, &lb, &send_elsize);
    MPI_Type_get_extent(recvtype, &lb, &recv_elsize);

#ifndef NO_ISEND_IRECV_IN_DOMAIN
    int n_requests;
    MPI_Request *requests = malloc(NTask * 2 * sizeof(MPI_Request));
    n_requests = 0;


    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;

        if(target >= NTask) continue;
        if(recvcnts[target] == 0) continue;
        MPI_Irecv(
                ((char*) recvbuf) + recv_elsize * rdispls[target], 
                recvcnts[target],
                recvtype, target, 101934, comm, &requests[n_requests++]);
    }

    MPI_Barrier(comm);
    /* not really necessary, but this will guarantee that all receives are
       posted before the sends, which helps the stability of MPI on
       bluegene, and perhaps some mpich1-clusters */
    /* Note 08/2016: Even on modern hardware this barrier leads to a slight speedup.
     * Probably because it allows the code to hit a fast path transfer.*/

    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;
        if(target >= NTask) continue;
        if(sendcnts[target] == 0) continue;
        MPI_Isend(((char*) sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target],
                sendtype, target, 101934, comm, &requests[n_requests++]);
    }

    MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
    free(requests);
#else
    for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
        int target = ThisTask ^ ngrp;

        if(target >= NTask) continue;
        if(sendcnts[target] == 0 && recvcnts[target] == 0) continue;
        MPI_Sendrecv(((char*)sendbuf) + send_elsize * sdispls[target], 
                sendcnts[target], sendtype, 
                target, 101934,
                ((char*)recvbuf) + recv_elsize * rdispls[target],
                recvcnts[target], recvtype, 
                target, 101934, 
                comm, MPI_STATUS_IGNORE);

    }
#endif
    /* ensure the collective-ness */
    MPI_Barrier(comm);

    return 0;
}

static void _mpsort_mpi_parse_env()
{
    static int _mpsort_env_parsed = 0;
    if(_mpsort_env_parsed) return;

    _mpsort_env_parsed = 1;
    if(getenv("MPSORT_ENABLE_SPARSE_ALLTOALLV"))
        mpsort_mpi_set_options(MPSORT_ENABLE_SPARSE_ALLTOALLV);
    if(getenv("MPSORT_DISABLE_IALLREDUCE"))
        mpsort_mpi_set_options(MPSORT_DISABLE_IALLREDUCE);
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
