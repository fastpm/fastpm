#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include "mpsort.h"

static double wtime() {
    struct timespec t1;
    clock_gettime(CLOCK_REALTIME, &t1);
    return (double)((t1.tv_sec+t1.tv_nsec*1e-9));
}

static void radix_int(const void * ptr, void * radix, void * arg) {
    *(int64_t*)radix = *(const int64_t*) ptr + INT64_MIN;
}
static int compar_int(const void * p1, const void * p2) {
    const int64_t * i1 = p1, *i2 = p2;
    return (*i1 > *i2) - (*i1 < *i2);
}

static int64_t
checksum(int64_t * data, size_t localsize, MPI_Comm comm)
{
    int64_t sum = 0;
    ptrdiff_t i;
    for(i = 0; i < localsize; i ++) {
        sum += data[i];
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_LONG, MPI_SUM, comm);
    return sum;
}
static void
generate(int64_t * data, size_t localsize, int bits, int seed)
{
    /* only keep bits of precision. */
    srandom(seed);

    ptrdiff_t i;
    unsigned shift = 64u - bits;
    for(i = 0; i < localsize; i ++) {
        uint64_t value = (int64_t) random() * (int64_t) random() * random() * random();
        data[i] = (signed) ((value << shift));
    }
}

static void
check_sorted(int64_t * data, size_t localsize, MPI_Comm comm)
{
    ptrdiff_t i;
    int ThisTask, NTask;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    int TAG = 0xbeef;

    ptrdiff_t dummy[1] = {INT64_MIN};

    for(i = 1; i < localsize; i ++) {
        if(data[i] < data[i - 1]) {
            fprintf(stderr, "Ordering of local array is broken. \n");
            MPI_Abort(comm, -1);
        }
    }

    if(NTask == 1) return;


    int64_t prev = -1;

    while(1) {
        if(ThisTask == 0) {
            int64_t * ptr = dummy;
            if(localsize > 0) {
                ptr = &data[localsize - 1];
            }
            MPI_Send(ptr, 1, MPI_LONG, ThisTask + 1, 0xbeef, comm);
            break;
        }
        if(ThisTask == NTask - 1) {
            MPI_Recv(&prev, 1, MPI_LONG,
                    ThisTask - 1, 0xbeef, comm, MPI_STATUS_IGNORE);
            break;
        }
        /* else */
        if(localsize == 0) {
            /* simply pass through whatever we get */
            MPI_Recv(&prev, 1, MPI_LONG, ThisTask - 1, 0xbeef, comm, MPI_STATUS_IGNORE);
            MPI_Send(&prev, 1, MPI_LONG, ThisTask + 1, 0xbeef, comm);
            break;
        }
        /* else */
        { 
            MPI_Sendrecv(
                    &data[localsize - 1], 1, MPI_LONG, 
                    ThisTask + 1, 0xbeef, 
                    &prev, 1, MPI_LONG,
                    ThisTask - 1, 0xbeef, comm, MPI_STATUS_IGNORE);
            break;
        }
    }

    if(ThisTask > 1) {
        if(localsize > 0) {
//                printf("ThisTask = %d prev = %d\n", ThisTask, prev);
            if(prev > data[0]) {
                fprintf(stderr, "global ordering fail\n");
                MPI_Abort(comm, -1);
            }
        }
    }
}
static void usage()
{
        printf("./main [number of items]\n");
}

int main(int argc, char * argv[]) {
    int i;
    MPI_Init(&argc, &argv);

    int ThisTask;
    int NTask;

    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    extern char * optarg;
    extern int optind;

    int opt;

    int staggered = 0;
    int bits;

    while(-1 != (opt = getopt(argc, argv, "IiSsGgb:"))) {
        switch(opt) {
            case 'b':
                bits = atoi(optarg);
                if(ThisTask == 0) {
                    printf("Bit width = %d\n", bits);
                }
                break;
            case 'i':
                mpsort_mpi_set_options(MPSORT_DISABLE_IALLREDUCE);
                if(ThisTask == 0) {
                    printf("DISABLE_IALLREDUCE\n");
                }
                break;
            case 'I':
                /* default */
                if(ThisTask == 0) {
                    printf("REQUIRE_IALLREDUCE\n");
                }
                break;
            case 'g':
                mpsort_mpi_set_options(MPSORT_DISABLE_GATHER_SORT);
                if(ThisTask == 0) {
                    printf("DISABLE_GATHER_SORT\n");
                }
                break;
            case 'G':
                mpsort_mpi_set_options(MPSORT_REQUIRE_GATHER_SORT);
                if(ThisTask == 0) {
                    printf("REQUIRE_GATHER_SORT\n");
                }
                break;
            case 's':
                /* Set some ranks to zero load */
                staggered = 1; 
                if(ThisTask == 0) {
                    printf("STAGGERED\n");
                }
                break;
            case 'S':
                /* Set some ranks to zero load */
                staggered = 0; 
                if(ThisTask == 0) {
                    printf("NOT STAGGERED\n");
                }
                break;
            default:
                usage();
                exit(-1);
        }
    }

    if(optind >= argc) {
        usage();
        exit(-1);
    }
    
    int64_t srcsize = atoi(argv[optind]);

    if(ThisTask == 0) {
        printf("NTask = %d\n", NTask);
        printf("src size = %ld\n", srcsize);
    }
    if(staggered && (ThisTask % 2 == 0)) srcsize = 0;


    int64_t csize;

    MPI_Allreduce(&srcsize, &csize, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    int64_t destsize = csize * (ThisTask + 1) /  NTask - csize * (ThisTask) / NTask;

    if(ThisTask == 0) {
        printf("dest size = %ld\n", destsize);
        printf("csize = %ld\n", csize);
    }
    int64_t * src = malloc(srcsize * sizeof(int64_t));
    int64_t * dest = malloc(destsize * sizeof(int64_t));

    int seed = 9999 * ThisTask;
    generate(src, srcsize, bits, seed);

    int64_t srcsum = checksum(src, srcsize, MPI_COMM_WORLD);

    {
        double start = MPI_Wtime();

        mpsort_mpi_newarray(src, srcsize,
                            dest, destsize,
                            sizeof(int64_t),
                            radix_int, sizeof(int64_t), NULL,
                            MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        double end = MPI_Wtime();

        int64_t destsum = checksum(dest, destsize, MPI_COMM_WORLD);

        if(destsum != srcsum) {
            printf("MPSort checksum is inconsistent.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        check_sorted(dest, destsize, MPI_COMM_WORLD);

        if(ThisTask == 0) {
            printf("MPSort total time: %g\n", end - start);
            mpsort_mpi_report_last_run();
        }
    }

    free(src);
    free(dest);

    MPI_Finalize();
    return 0;
}
