#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include "mpsort.h"

static void radix_int(const void * ptr, void * radix, void * arg) {
    *(int*)radix = *(const int*) ptr;
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
    int * mydata = malloc(mysize * sizeof(int));

    int64_t mysum = 0;
    int64_t truesum = 0, realsum = 0;

    for(i = 0; i < mysize; i ++) {
        mydata[i] = random() % 10000;
        mysum += mydata[i];
    }

    MPI_Allreduce(&mysum, &truesum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    mpsort_mpi(mydata, mysize, sizeof(int),
            radix_int, sizeof(int),
            NULL, MPI_COMM_WORLD);

    if(ThisTask == 0)  {
        mpsort_mpi_report_last_run();
    }

    mysum = 0;
    for(i = 0; i < mysize; i ++) {
        mysum += mydata[i];
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
        int prev = -1;
        if(ThisTask == 0) {
            if(mysize == 0) {
                MPI_Send(&prev, 1, MPI_INT, 
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            } else {
                MPI_Send(&mydata[mysize - 1], 1, MPI_INT, 
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            }
        } else
        if(ThisTask == NTask - 1) {
            MPI_Recv(&prev, 1, MPI_INT,
                    ThisTask - 1, 0xbeef, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            if(mysize == 0) {
                MPI_Recv(&prev, 1, MPI_INT,
                        ThisTask - 1, 0xbeef, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&prev, 1, MPI_INT,
                        ThisTask + 1, 0xbeef, MPI_COMM_WORLD);
            } else {
                MPI_Sendrecv(
                        &mydata[mysize - 1], 1, MPI_INT, 
                        ThisTask + 1, 0xbeef, 
                        &prev, 1, MPI_INT,
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
