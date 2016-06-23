#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include "mpsort.h"

#include "parallel_sort.c"

static double wtime() {
    struct timespec t1;
    clock_gettime(CLOCK_REALTIME, &t1);
    return (double)((t1.tv_sec+t1.tv_nsec*1e-9));
}

static void radix_int(const void * ptr, void * radix, void * arg) {
    *(int64_t*)radix = *(const int64_t*) ptr;
}
static int compar_int(const void * p1, const void * p2) {
    const int64_t * i1 = p1, *i2 = p2;
    return (*i1 > *i2) - (*i1 < *i2);
}

int main(int argc, char * argv[]) {
    int i;
    MPI_Init(&argc, &argv);

    int ThisTask;
    int NTask;

    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    srand(9999 * ThisTask);

    if(argc != 2) {
        printf("./main [number of items]\n");
        return 1;
    }

    size_t mysize = atoi(argv[1]);
    int64_t * mydata = malloc(mysize * sizeof(int64_t));
    int64_t * mydata2 = malloc(mysize * sizeof(int64_t));

    int64_t mysum = 0;
    int64_t truesum = 0, realsum = 0;

    if(ThisTask == 2) mysize = 0;

    for(i = 0; i < mysize; i ++) {
        uint64_t data = (int64_t) random() * (int64_t) random() * random() * random();
        //data = 0 * ThisTask * (int64_t) mysize + i / 10;
        mydata[i] = data & 0xffffffffffffff;
        mydata2[i] = mydata[i];
        mysum += mydata[i];
    }

    MPI_Allreduce(&mysum, &truesum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    {
        double start = MPI_Wtime();
        radix_sort_mpi(mydata, mysize, sizeof(int64_t),
                radix_int, sizeof(int64_t),
                NULL, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        double end = MPI_Wtime();
        if(ThisTask == 0)
            printf("radix sort: %g\n", end - start);
    }
    if(ThisTask == 0)  {
        radix_sort_mpi_report_last_run();
    }

    {
        double start = MPI_Wtime();
        parallel_sort(mydata2, mysize, sizeof(int64_t),
                compar_int);
        MPI_Barrier(MPI_COMM_WORLD);
        double end = MPI_Wtime();
        if(ThisTask == 0)
        printf("parallel sort: %g\n", end - start);
    }
    mysum = 0;
    for(i = 0; i < mysize; i ++) {
        mysum += mydata[i];
        if(mydata[i] != mydata2[i]) {
            fprintf(stderr, "sorting error\n");
            abort();
        }
    }

    MPI_Allreduce(&mysum, &realsum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(realsum != truesum) {
        fprintf(stderr, "checksum fail\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for(i = 1; i < mysize; i ++) {
        if(mydata[i] < mydata[i - 1]) {
            fprintf(stderr, "local ordering fail\n");
        }
    }
    if(NTask > 1) {
        int64_t prev = -1;
        if(ThisTask == 0) {
            if(mysize == 0) {
                MPI_Send(&prev, 1, MPI_LONG_LONG, 
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            } else {
                MPI_Send(&mydata[mysize - 1], 1, MPI_LONG_LONG, 
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            }
        } else
        if(ThisTask == NTask - 1) {
            MPI_Recv(&prev, 1, MPI_LONG_LONG,
                    ThisTask - 1, 0xbeef, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            if(mysize == 0) {
                MPI_Recv(&prev, 1, MPI_LONG_LONG,
                        ThisTask - 1, 0xbeef, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&prev, 1, MPI_LONG_LONG,
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            } else {
                MPI_Sendrecv(
                        &mydata[mysize - 1], 1, MPI_LONG_LONG, 
                        ThisTask + 1, 0xbeef, 
                        &prev, 1, MPI_LONG_LONG,
                        ThisTask - 1, 0xbeef, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if(ThisTask > 1) {
            if(mysize > 0) {
//                printf("ThisTask = %d prev = %d\n", ThisTask, prev);
                if(prev > mydata[0]) {
                    fprintf(stderr, "global ordering fail\n");
                    abort();
                }
            }
        }
    }
    free(mydata);
    MPI_Finalize();
    return 0;
}
