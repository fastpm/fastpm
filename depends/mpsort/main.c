#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "mpsort.h"

static double wtime() {
    struct timespec t1;
    clock_gettime(CLOCK_REALTIME, &t1);
    return (double)((t1.tv_sec+t1.tv_nsec*1e-9));
}

static void radix_int(const void * ptr, void * radix, void * arg) {
    *(int*)radix = *(const int*) ptr;
}
static int compar_int(const void * p1, const void * p2) {
    const unsigned int * i1 = p1, *i2 = p2;
    return (*i1 > *i2) - (*i1 < *i2);
}

int main(int argc, char * argv[]) {
    int i;
    srand(9999);
    if(argc != 2) {
        printf("./main [number of items]\n");
        return 1;
    }
    int NUMITEMS = atoi(argv[1]);
    int * data1 = malloc(sizeof(int) * NUMITEMS);
    int * data2 = malloc(sizeof(int) * NUMITEMS);
    int * data3 = malloc(sizeof(int) * NUMITEMS);
    for(i = 0; i < NUMITEMS; i ++) {
        data1[i] = random() % 10000;
        data2[i] = data1[i];
        data3[i] = data1[i];
    }

    {
        double t0 = wtime();
        radix_sort_omp(data2, NUMITEMS, sizeof(int),
                radix_int, sizeof(int),
                NULL);
        double t1 = wtime();
        printf("time spent omp: %g\n", t1 - t0);
    }

    {
        double t0 = wtime();
        qsort(data3, NUMITEMS, sizeof(int), compar_int);
        double t1 = wtime();
        printf("time spent qsort: %g\n", t1 - t0);
    }
    {
        double t0 = wtime();
        radix_sort(data1, NUMITEMS, sizeof(int),
                radix_int, sizeof(int),
                NULL);
        printf("max is %u\n", data1[NUMITEMS - 1]);
        double t1 = wtime();
        printf("time spent radix: %g\n", t1 - t0);
    }
    for(i = 0; i < NUMITEMS; i ++) {
        if(data1[i] != data2[i]) abort();
    }
    free(data1);
    free(data2);
    free(data3);
    return 0;
}
