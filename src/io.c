#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
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

    FastPMSolverBase * base = &fastpm->base;

    int ThisTask = base->ThisTask;
    int NTask = base->NTask;
    MPI_Comm comm = base->comm;

    if(Nwriters == 0 || Nwriters > NTask) Nwriters = NTask;

    int Nfile = NTask / 8;
    if (Nfile == 0) Nfile = 1;
    int64_t size = p->np;
    int64_t offsets[NTask];

    double H0 = 100.;
    double RSD = 1.0 / (H0 * p->a_x * HubbleEa(p->a_x, CP(fastpm)));

    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);

    int i;
    BigFile bf;
    big_file_mpi_create(&bf, filebase, comm);
    {
        BigBlock bb;
        big_file_mpi_create_block(&bf, &bb, "header", "i8", 0, 1, 0, comm);
        double ScalingFactor = p->a_x;
        double OmegaM = fastpm->omega_m;
        double BoxSize = fastpm->boxsize;
        uint64_t NC = fastpm->nc;
        double rho_crit = 27.7455;
        double M0 = OmegaM * rho_crit * (BoxSize / NC) * (BoxSize / NC) * (BoxSize / NC);

        big_block_set_attr(&bb, "BoxSize", &BoxSize, "f8", 1);
        big_block_set_attr(&bb, "ScalingFactor", &ScalingFactor, "f8", 1);
        big_block_set_attr(&bb, "RSDFactor", &RSD, "f8", 1);
        big_block_set_attr(&bb, "OmegaM", &OmegaM, "f8", 1);
        big_block_set_attr(&bb, "NC", &NC, "i8", 1);
        big_block_set_attr(&bb, "M0", &M0, "f8", 1);
        big_block_set_attr(&bb, "ParamFile", parameters, "S1", strlen(parameters) + 1);
        big_block_mpi_close(&bb, comm);
    }
    struct {
        char * name;
        void * base;
        char * dtype;
        int nmemb;
        char * dtype_out;
    } * bdesc, BLOCKS[] = {
        {"Position", p->x, "f8", 3, "f4"},
        {"Velocity", p->v, "f4", 3, "f4"},
        {"ID", p->id, "i8", 1, "i8"},
        {NULL, },
    };

    for(bdesc = BLOCKS; bdesc->name; bdesc ++) {
        fastpm_info("Writing block %s\n", bdesc->name);
        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        big_file_mpi_create_block(&bf, &bb, bdesc->name, bdesc->dtype_out, bdesc->nmemb,
                    Nfile, size, comm);

        big_array_init(&array, bdesc->base, bdesc->dtype, 2, (size_t[]) {p->np, bdesc->nmemb}, NULL);
        big_block_seek(&bb, &ptr, 0);
        big_block_mpi_write(&bb, &ptr, &array, Nwriters, comm);
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

