#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include "mpsort.h"

void radix(const void * ptr, void * radix, void * arg) {
    uint64_t * u = (uint64_t *) radix;
    memset(u, 0, 16);
    memcpy(&u[0], ((const char*) ptr) + 16, 8);
    memcpy(&u[1], ((const char*) ptr) + 0, 4);
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    int ThisTask;
    int NTask;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    char fnbuf[1024];
    sprintf(fnbuf, "test-data/mpsort_mpi.b40.r16.0-%05d-of-%05d", ThisTask, NTask);

    FILE * fp = fopen(fnbuf, "r");
    fseek(fp, 0, SEEK_END);
    size_t sz = ftell(fp);

    void * base = malloc(sz);
    fseek(fp, 0, SEEK_SET);
    fread(base, 40, sz, fp);
    fclose(fp);

    mpsort_mpi_set_options(MPSORT_DISABLE_GATHER_SORT);
    mpsort_mpi(base, sz / 40, 40, radix, 16, NULL, comm);

    int i;
    int r;
    for(r = 0; r < ThisTask; r ++) {
        MPI_Barrier(comm);
    }
    for(i = 0; i < sz / 40; i ++) {
        uint64_t r[2];
        radix(((char*) base) + i * 40, r, NULL);
        printf("%d : %ld %ld \n", ThisTask, r[0], r[1]);
    }
    for(; r < NTask; r ++) {
        MPI_Barrier(comm);
    }
    free(base);
}
