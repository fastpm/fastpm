#include <mpi.h>
#include <stdio.h>
#include <bigfile.h>
#include <bigfile-mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>
static void 
cumsum(int64_t offsets[], int N) 
{
    int i;
    int64_t tmp[N];

    tmp[0] = 0;
    for(i = 1; i < N; i ++) {
        tmp[i] = tmp[i - 1] + offsets[i - 1];
    }
    for(i = 0; i < N; i ++) {
        offsets[i] = tmp[i];
    }
}

int 
write_snapshot(FastPM * fastpm, PMStore * p, char * filebase) 
{
    MPI_Comm comm = fastpm->comm;
    int NTask;
    int ThisTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    int Nfile = NTask / 32;
    if (Nfile == 0) Nfile = 1;

    int64_t size = p->np;
    int64_t offsets[NTask];

    MPI_Allgather(&size, 1, MPI_LONG, offsets, 1, MPI_LONG, comm);
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);
    cumsum(offsets, NTask);
 
    BigFile bf = {0};
    big_file_mpi_create(&bf, filebase, comm);

    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "Position", "f8", 3, Nfile, size, comm);
        big_array_init(&array, p->x, "f8", 2, (size_t[]) {p->np, 3}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        big_block_write(&bb, &ptr, &array);
        big_block_mpi_close(&bb, comm);
    }
    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "Velocity", "f4", 3, Nfile, size, comm);
        big_array_init(&array, p->v, "f4", 2, (size_t[]) {p->np, 3}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        big_block_write(&bb, &ptr, &array);
        big_block_mpi_close(&bb, comm);
    }
    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "ID", "i8", 1, Nfile, size, comm);
        big_array_init(&array, p->id, "i8", 2, (size_t[]) {p->np, 1}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        big_block_write(&bb, &ptr, &array);
        big_block_mpi_close(&bb, comm);
    }
    big_file_mpi_close(&bf, comm);
    return 0;
}

int 
read_snapshot(FastPM * fastpm, PMStore * p, char * filebase) 
{
    return 0;
}

