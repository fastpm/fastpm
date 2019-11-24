#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "internal.h"

#include "internal-parallel.h"
#include "mpsort.h"

struct crompstruct {
    struct piter pi;
    unsigned char * P;
    unsigned char * Pmax;
    unsigned char * Pmin;
    ptrdiff_t * C; /* expected counts */
    ptrdiff_t * CLT; /* counts of less than P */
    ptrdiff_t * CLE; /* counts of less than or equal to P */
    ptrdiff_t * GL_CLT; /* counts of less than P */
    ptrdiff_t * GL_CLE; /* counts of less than or equal to P */
    ptrdiff_t * GL_C; /* counts of less than or equal to P */
};


/* OPENMP version of radix sort;
 * this is truely parallel;
 * but it is usually slower than
 * simple radix sort if the number of processor is small.
 *
 * some benchmarks on Coma at CMU shows best performance is at 16
 * CPUs; still faster than serial version with 8 CPUs.
 * comparable with qsort (non-radix) at 8 CPUs.
 *
 * the coding is more of a prototype of the MPI radix sort;
 * it is thus very poorly written in the standards of an OPENMP program;
 * */

static void _setup_mpsort_omp(struct crompstruct * o, struct crstruct * d);
static void _cleanup_mpsort_omp(struct crompstruct * o, struct crstruct * d);

static void mpsort_omp_single(void * base, size_t nmemb,
        struct crstruct * d, struct crompstruct * o);

void mpsort_omp(void * base, size_t nmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg) {
    if(nmemb == 0) return;

    struct crstruct d;
    struct crompstruct o;
    _setup_radix_sort(&d, base, nmemb, size, radix, rsize, arg);

    _setup_mpsort_omp(&o, &d);

    /*
     * first solve for P such that CLT[i] < C <= CLE[i]
     *
     * Then calculate a communication layout.
     *
     * Then alltoall.
     * Then local sort again
     * */

#pragma omp parallel
    {
        mpsort_omp_single (base, nmemb, &d, &o);
    }

    _cleanup_mpsort_omp(&o, &d);
}

static void _setup_mpsort_omp(struct crompstruct * o, struct crstruct * d) {
    int NTaskMax = omp_get_max_threads();

    o->P = calloc(NTaskMax, d->rsize);

    int NTaskMax1 = NTaskMax + 1;
    size_t NTaskMax12 = (NTaskMax + 1) * (NTaskMax + 1);

    /* following variables are used in counting index by ThisTask + 1 */
    o->GL_CLT = calloc(NTaskMax12, sizeof(ptrdiff_t));
    o->GL_CLE = calloc(NTaskMax12, sizeof(ptrdiff_t));
    o->GL_C = calloc(NTaskMax12, sizeof(ptrdiff_t));

    o->C = calloc(NTaskMax1, sizeof(ptrdiff_t)); /* expected counts */
    o->CLT = calloc(NTaskMax1, sizeof(ptrdiff_t)); /* counts of less than P */
    o->CLE = calloc(NTaskMax1, sizeof(ptrdiff_t)); /* counts of less than or equal to P */

    o->Pmin = calloc(1, d->rsize);
    o->Pmax = calloc(1, d->rsize);
    memset(o->Pmin, -1, d->rsize);
    memset(o->Pmax, 0, d->rsize);
}
static void _cleanup_mpsort_omp(struct crompstruct * o, struct crstruct * d) {
    free(o->CLE);
    free(o->CLT);
    free(o->C);
    free(o->GL_C);
    free(o->GL_CLE);
    free(o->GL_CLT);
    free(o->Pmin);
    free(o->Pmax);
    free(o->P);
}

static void _reduce_sum(ptrdiff_t * send, ptrdiff_t * recv, size_t count) {
    ptrdiff_t i;
#pragma omp single
    for(i = 0; i < count; i ++) {
        recv[i] = 0;
    }
#pragma omp critical
    for(i = 0; i < count; i ++) {
        recv[i] += send[i];
    }
#pragma omp barrier
}
static void _gather(void * sendbuf, int sendcount1, void * recvbuf, size_t itemsize) {
    int ThisTask = omp_get_thread_num();
    memcpy((char*) recvbuf + ThisTask * sendcount1 * itemsize,
            sendbuf, sendcount1 * itemsize);
#pragma omp barrier
}

/*
 * solve for the communication layout based on
 *
 * C: the desired number of items per task
 * GL_CLT[t,i+1]: the offset of lt P[i] in task t
 * GL_CLE[t,i+1]: the offset of le P[i] in task t
 *
 * the result is saved in
 *
 * GL_C[t, i]: the offset of sending to task i in task t.
 *
 * this routine requires GL_ to scale with NTask * NTask;
 * won't work with 1,000 + ranks.
 * */
static void _solve_for_layout (
        int NTask,
        ptrdiff_t * C,
        ptrdiff_t * GL_CLT,
        ptrdiff_t * GL_CLE,
        ptrdiff_t * GL_C) {
    int NTask1 = NTask + 1;
    int i, j;
    /* first assume we just send according to GL_CLT */
    for(i = 0; i < NTask + 1; i ++) {
        for(j = 0; j < NTask; j ++) {
            GL_C[j * NTask1 + i] = GL_CLT[j * NTask1 + i];
        }
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

    for(i = 0; i < NTask; i ++) {
        ptrdiff_t sure = 0;

        /* how many will I surely receive? */
        for(j = 0; j < NTask; j ++) {
            ptrdiff_t sendcount = GL_C[j * NTask1 + i + 1] - GL_C[j * NTask1 + i];
            sure += sendcount;
        }
        /* let's see if we have enough */
        ptrdiff_t deficit = C[i + 1] - C[i] - sure;

        for(j = 0; j < NTask; j ++) {
            /* deficit solved */
            if(deficit == 0) break;
            if(deficit < 0) {
                fprintf(stderr, "serious bug: more items than there should be: deficit=%ld\n", deficit);
                abort();
            }
            /* how much task j can supply ? */
            ptrdiff_t supply = GL_CLE[j * NTask1 + i + 1] - GL_C[j * NTask1 + i + 1];
            if(supply < 0) {
                fprintf(stderr, "serious bug: less items than there should be: supply =%ld\n", supply);
                abort();
            }
            if(supply <= deficit) {
                GL_C[j * NTask1 + i + 1] += supply;
                deficit -= supply;
            } else {
                GL_C[j * NTask1 + i + 1] += deficit;
                deficit = 0;
            }
        }
    }

#if 0
    for(i = 0; i < NTask; i ++) {
        for(j = 0; j < NTask + 1; j ++) {
            printf("%d %d %d, ",
                    GL_CLT[i * NTask1 + j],
                    GL_C[i * NTask1 + j],
                    GL_CLE[i * NTask1 + j]);
        }
        printf("\n");
    }
#endif

}

static void mpsort_omp_single(void * base, size_t nmemb,
        struct crstruct * d, struct crompstruct * o) {
    int NTask = omp_get_num_threads();
    int ThisTask = omp_get_thread_num();

    ptrdiff_t myCLT[NTask + 1]; /* counts of less than P */
    ptrdiff_t myCLE[NTask + 1]; /* counts of less than or equal to P */

    int i;

#pragma omp single
    {
        o->C[0] = 0;
        for(i = 0; i < NTask; i ++) {
        /* how many items are desired per thread */
            o->C[i + 1] = nmemb * (i + 1) / NTask;
        }
    }

    double t0 = omp_get_wtime();

    /* distribute the array evenly */
    char * mybase = (char*) base + nmemb * ThisTask / NTask * d->size;
    size_t mynmemb = nmemb * (ThisTask + 1)/ NTask - nmemb * (ThisTask) / NTask;


    /* and sort the local array */
    radix_sort(mybase, mynmemb, d->size, d->radix, d->rsize, d->arg);


    /* find the max radix and min radix of all */
    if(mynmemb > 0) {
        unsigned char myPmax[d->rsize];
        unsigned char myPmin[d->rsize];
        d->radix(mybase + (mynmemb - 1) * d->size, myPmax, d->arg);
        d->radix(mybase, myPmin, d->arg);
#pragma omp critical
        {
            if(d->compar(myPmax, o->Pmax, d->rsize) > 0) {
                memcpy(o->Pmax, myPmax, d->rsize);
            }
            if(d->compar(myPmin, o->Pmin, d->rsize) < 0) {
                memcpy(o->Pmin, myPmin, d->rsize);
            }
        }
    }
#pragma omp barrier

    double t1 = omp_get_wtime();
    printf("Initial sort took %g\n", t1 - t0);

    /* now do the radix counting iterations */

#pragma omp single
    piter_init(&o->pi, o->Pmin, o->Pmax, NTask - 1, d);

    int iter = 0;

    int done = 0;

    while(!done) {
        iter ++;
#pragma omp barrier
#pragma omp single
        piter_bisect(&o->pi, o->P);

        _histogram(o->P, NTask - 1, mybase, mynmemb, myCLT, myCLE, d);

        _reduce_sum(myCLT, o->CLT, NTask + 1);
        _reduce_sum(myCLE, o->CLE, NTask + 1);

#pragma omp single
        piter_accept(&o->pi, o->P, o->C, o->CLT, o->CLE);

        done = piter_all_done(&o->pi);
    }
#pragma omp barrier

#pragma omp single
    piter_destroy(&o->pi);

    double t2 = omp_get_wtime();
    printf("counting took %g\n", t2 - t1);

#if 0
#pragma omp single
    {
        printf("AfterIteration: split , CLT, C, CLE\n");
        int i;
        for(i = 0; i < NTask + 1; i ++) {
            printf("%d %ld %ld %ld\n", i, o->CLT[i], o->C[i], o->CLE[i]);
        }
    }
#endif

    _histogram(o->P, NTask - 1, mybase, mynmemb, myCLT, myCLE, d);

    /* gather to all (used only by single */
    _gather(myCLT, NTask + 1, o->GL_CLT, sizeof(ptrdiff_t));
    _gather(myCLE, NTask + 1, o->GL_CLE, sizeof(ptrdiff_t));

    /* indexing is
     * GL_C[Sender * NTask1 + Recver] */
    /* here we know:
     * o->GL_CLT, o->GL_CLE
     * The first NTask - 1 items in CLT and CLE gives the bounds of
     * split points  ( CLT < split <= CLE)
     * */

    /* find split points in O->GL_C*/
#pragma omp single
    _solve_for_layout(NTask, o->C, o->GL_CLT, o->GL_CLE, o->GL_C);

    double t3 = omp_get_wtime();
    printf("find split took %g\n", t3 - t2);

    /* exchange data */
    /* */
    char * buffer = malloc(d->size * mynmemb);
    int NTask1 = NTask + 1;

#if 0
#pragma omp critical
    {
        printf("%d contains %d items ", ThisTask, mynmemb);
        int k;
        int * ptr = mybase;
        for(k = 0; k < mynmemb; k ++) {
            printf("%d ", ((int *) ptr)[k]);
        }
        printf("\n");
    }
#pragma omp barrier
#endif

    /* can't do with an alltoall because it is difficult to sync
     * the sendbuf pointers with omp */
    char * recv = buffer;
    for(i = 0; i < NTask; i ++) {
        char * sendbuf = (char*) base + d->size * (i * nmemb / NTask);
        size_t size = (o->GL_C[i * NTask1 + ThisTask + 1]
                - o->GL_C[i * NTask1 + ThisTask]) * d->size;
        char * ptr = &sendbuf[d->size * o->GL_C[i * NTask1 + ThisTask]];

        memcpy(recv, ptr, size);
        recv += size;
    }

#pragma omp barrier
    memcpy(mybase, buffer, mynmemb * d->size);
    free(buffer);

#if 0
#pragma omp critical
    {
        printf("%d after exchange %d items ", ThisTask, mynmemb);
        int k;
        int * ptr = mybase;
        for(k = 0; k < mynmemb; k ++) {
            printf("%d ", ((int *) ptr)[k]);
        }
        printf("\n");
    }
#pragma omp barrier
#endif

    double t4 = omp_get_wtime();
    printf("exchange took %g\n", t4 - t3);

    radix_sort(mybase, mynmemb, d->size, d->radix, d->rsize, d->arg);

    double t5 = omp_get_wtime();
    printf("final sort %g\n", t5 - t4);

#pragma omp barrier
}

