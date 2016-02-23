#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <bigfile.h>
#include <bigfile-mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
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

static Cosmology CP(FastPM * fastpm) {
    Cosmology c = {
        .OmegaM = fastpm->omega_m,
        .OmegaLambda = 1 - fastpm->omega_m,
    };
    return c;
}

int 
write_snapshot(FastPM * fastpm, PMStore * p, char * filebase, char * parameters, int Nwriters) 
{
    MPI_Comm comm = fastpm->comm;
    int NTask;
    int ThisTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    if(Nwriters == 0 || Nwriters > NTask) Nwriters = NTask;

    int Nfile = NTask / 8;
    if (Nfile == 0) Nfile = 1;
    int64_t size = p->np;
    int64_t offsets[NTask];

    double vfac = 100. / p->a_x;
    double RSD = p->a_x / Qfactor(p->a_x, CP(fastpm)) / vfac;

    MPI_Allgather(&size, 1, MPI_LONG, offsets, 1, MPI_LONG, comm);
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);
    cumsum(offsets, NTask);

    int i;
    BigFile bf = {0};
    big_file_mpi_create(&bf, filebase, comm);
    {
        BigBlock bb = {0};
        big_file_mpi_create_block(&bf, &bb, "header", "i8", 0, 1, 0, comm);
        double ScalingFactor = p->a_x;
        double OmegaM = fastpm->omega_m;
        double BoxSize = fastpm->boxsize;
        uint64_t NC = fastpm->nc;
        big_block_set_attr(&bb, "BoxSize", &BoxSize, "f8", 1);
        big_block_set_attr(&bb, "ScalingFactor", &ScalingFactor, "f8", 1);
        big_block_set_attr(&bb, "RSDFactor", &RSD, "f8", 1);
        big_block_set_attr(&bb, "OmegaM", &OmegaM, "f8", 1);
        big_block_set_attr(&bb, "NC", &NC, "i8", 1);
        big_block_set_attr(&bb, "ParamFile", parameters, "S1", strlen(parameters) + 1);
        big_block_mpi_close(&bb, comm);
    }
    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "Position", "f4", 3, Nfile, size, comm);
        big_array_init(&array, p->x, "f8", 2, (size_t[]) {p->np, 3}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        for(i = 0; i < Nwriters; i ++) {
            MPI_Barrier(comm);
            if(ThisTask % Nwriters != i) continue;
            big_block_write(&bb, &ptr, &array);
        }
        big_block_mpi_close(&bb, comm);
    }
    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "Velocity", "f4", 3, Nfile, size, comm);
        big_array_init(&array, p->v, "f4", 2, (size_t[]) {p->np, 3}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        for(i = 0; i < Nwriters; i ++) {
            MPI_Barrier(comm);
            if(ThisTask % Nwriters != i) continue;
            big_block_write(&bb, &ptr, &array);
        }
        big_block_mpi_close(&bb, comm);
    }
    {
        BigBlock bb = {0};
        BigArray array = {0};
        BigBlockPtr ptr = {0};
        big_file_mpi_create_block(&bf, &bb, "ID", "i8", 1, Nfile, size, comm);
        big_array_init(&array, p->id, "i8", 2, (size_t[]) {p->np, 1}, NULL);
        big_block_seek(&bb, &ptr, offsets[ThisTask]);
        for(i = 0; i < Nwriters; i ++) {
            MPI_Barrier(comm);
            if(ThisTask % Nwriters != i) continue;
            big_block_write(&bb, &ptr, &array);
        }
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

