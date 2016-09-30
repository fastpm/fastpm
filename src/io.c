#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bigfile.h>
#include <bigfile-mpi.h>
#include <mpsort.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>

/*
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
*/
int 
write_snapshot(FastPMSolver * fastpm, FastPMStore * p, char * filebase, char * parameters, int Nwriters) 
{

    int NTask = fastpm->NTask;
    MPI_Comm comm = fastpm->comm;

    if(Nwriters == 0 || Nwriters > NTask) Nwriters = NTask;

    int Nfile = NTask / 8;
    if (Nfile == 0) Nfile = 1;
    int64_t size = p->np;

    double H0 = 100.;
    /* Conversion from peculiar velocity to RSD,
     * http://mwhite.berkeley.edu/Talks/SantaFe12_RSD.pdf */
    double RSD = 1.0 / (H0 * p->a_x * HubbleEa(p->a_x, fastpm->cosmology));

    fastpm_info("RSD factor %e\n", RSD);

    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);

    BigFile bf;
    big_file_mpi_create(&bf, filebase, comm);
    {
        BigBlock bb;
        if(0 != big_file_mpi_create_block(&bf, &bb, ".", "i8", 0, 1, 0, comm)) {
            fastpm_raise(-1, "Failed to create the attributes\n");
        }
        double ScalingFactor = p->a_x;
        double OmegaM = fastpm->cosmology->OmegaM;
        double BoxSize = fastpm->config->boxsize;
        uint64_t NC = fastpm->config->nc;
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
        void * fastpm;
        char * dtype;
        int nmemb;
        char * dtype_out;
    } * bdesc, BLOCKS[] = {
        {"Position", p->x, "f8", 3, "f4"},
        {"Velocity", p->v, "f4", 3, "f4"},
        {"ID", p->id, "i8", 1, "i8"},
        {"Aemit", p->aemit, "f4", 1, "f4"},
        {NULL, },
    };

    for(bdesc = BLOCKS; bdesc->name; bdesc ++) {
        if(bdesc->fastpm == NULL) continue;

        fastpm_info("Writing block %s\n", bdesc->name);
        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        if(0 != big_file_mpi_create_block(&bf, &bb, bdesc->name, bdesc->dtype_out, bdesc->nmemb,
                    Nfile, size, comm)) {
            fastpm_raise(-1, "Failed to create the block\n");
        }

        big_array_init(&array, bdesc->fastpm, bdesc->dtype, 2, (size_t[]) {p->np, bdesc->nmemb}, NULL);
        big_block_seek(&bb, &ptr, 0);
        big_block_mpi_write(&bb, &ptr, &array, Nwriters, comm);
        big_block_mpi_close(&bb, comm);
    }

    big_file_mpi_close(&bf, comm);
    return 0;
}

int 
read_snapshot(FastPMSolver * fastpm, FastPMStore * p, char * filebase)
{
    return 0;
}

struct BufType {
    uint64_t ind;
    uint64_t ind2;
    float value[2];
};

void _radix(const void * ptr, void * radix, void * arg)
{
    const struct BufType * buf = ptr;
    uint64_t * key = radix;
    *key = buf->ind;
}
void _radix2(const void * ptr, void * radix, void * arg)
{
    const struct BufType * buf = ptr;
    uint64_t * key = radix;
    *key = buf->ind2;
}

int
write_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters)
{
    MPI_Comm comm = pm_comm(pm);
    BigFile bf;
    PMKIter kiter;
    if(Nwriters == 0) {
        MPI_Comm_size(comm, &Nwriters);
    }
    int ThisTask;
    int NTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    struct BufType * buf = malloc(sizeof(struct BufType) * pm_allocsize(pm) / 2);

    int Nmesh = pm_nmesh(pm)[0];
    double BoxSize = pm_boxsize(pm)[0];
    int64_t strides[3] = {Nmesh * (Nmesh / 2 + 1), Nmesh / 2 + 1, 1};
    int64_t shape[3] = {Nmesh , Nmesh, Nmesh / 2 + 1};

    size_t localsize = 0;
    for(pm_kiter_init(pm, &kiter);
       !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        uint64_t iabs = kiter.iabs[0] * strides[0] + kiter.iabs[1] * strides[1] + kiter.iabs[2] * strides[2];
        buf[localsize].value[0] = data[kiter.ind];
        buf[localsize].value[1] = data[kiter.ind + 1];
        buf[localsize].ind2 = iabs;
        buf[localsize].ind = ThisTask * ((size_t) Nmesh) * Nmesh * Nmesh + localsize;
        localsize ++;
    }

    /* sort by k */
    mpsort_mpi(buf, localsize, sizeof(buf[0]), _radix2, 8, NULL, comm);

    size_t size = localsize;
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);


    int Nfile = NTask / 8;
    if (Nfile == 0) Nfile = 1;

    big_file_mpi_create(&bf, filename, comm);
    {
        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        big_file_mpi_create_block(&bf, &bb, blockname, "c8", 1, Nfile, size, comm);
        big_array_init(&array, &buf[0].value[0], "c8", 1, (size_t[]) {localsize, 1},
                (ptrdiff_t[]) { sizeof(buf[0]), sizeof(buf[0])});
        big_block_seek(&bb, &ptr, 0);
        big_block_mpi_write(&bb, &ptr, &array, Nwriters, comm);

        big_block_set_attr(&bb, "ndarray.ndim", (int[]){3,}, "i4", 1);
        big_block_set_attr(&bb, "ndarray.strides", strides, "i8", 3);
        big_block_set_attr(&bb, "ndarray.shape", shape, "i8", 3);
        big_block_set_attr(&bb, "Nmesh", &Nmesh, "i4", 1);
        big_block_set_attr(&bb, "BoxSize", &BoxSize, "f8", 1);
        big_block_mpi_close(&bb, comm);
    }

    free(buf);
    big_file_mpi_close(&bf, comm);
    return 0;
}

int
read_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters)
{
    MPI_Comm comm = pm_comm(pm);
    BigFile bf;
    PMKIter kiter;
    if(Nwriters == 0) {
        MPI_Comm_size(comm, &Nwriters);
    }

    struct BufType * buf = malloc(sizeof(struct BufType) * pm_allocsize(pm) / 2);

    int Nmesh = pm_nmesh(pm)[0];
    ptrdiff_t strides[3] = {Nmesh * (Nmesh / 2 + 1), Nmesh / 2 + 1, 1};
    int64_t shape[3] = {Nmesh , Nmesh, Nmesh / 2 + 1};

    int ThisTask;
    int NTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    uint64_t localsize = 0;
    for(pm_kiter_init(pm, &kiter);
       !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        uint64_t iabs = kiter.iabs[0] * strides[0] + kiter.iabs[1] * strides[1] + kiter.iabs[2] * strides[2];
        buf[localsize].value[0] = 0;
        buf[localsize].value[1] = 0;
        buf[localsize].ind2 = iabs;
        buf[localsize].ind = ThisTask * ((size_t) Nmesh) * Nmesh * Nmesh + localsize;
        localsize ++;
    }
    /* sort by ind2, such that ind is the original linear location in 2d decomposition */
    mpsort_mpi(buf, localsize, sizeof(buf[0]), _radix2, 8, NULL, comm);

    int Nfile = NTask / 8;
    if (Nfile == 0) Nfile = 1;

    big_file_mpi_open(&bf, filename, comm);
    {
        int64_t istrides[3];
        int64_t ishape[3];

        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        big_file_mpi_open_block(&bf, &bb, blockname, comm);

        big_block_get_attr(&bb, "ndarray.strides", istrides, "i8", 3);
        big_block_get_attr(&bb, "ndarray.shape", ishape, "i8", 3);

        int d;
        for(d = 0; d < 3; d++) {
            if(ishape[d] != shape[d]) {
                fastpm_raise(-1, "Shape of complex field mismatch. Expecting (%ld %ld %ld), file has (%ld %ld %ld)\n",
                    shape[0], shape[1], shape[2], ishape[0], ishape[1], ishape[2]);
            }
        }
        /* FIXME: assert strides is consistent */
        big_array_init(&array, &buf[0].value[0], "c8", 1, (size_t[]) {localsize, 1}, (ptrdiff_t[]) { sizeof(buf[0]), sizeof(buf[0])});
        big_block_seek(&bb, &ptr, 0);
        big_block_mpi_read(&bb, &ptr, &array, Nwriters, comm);

        big_block_mpi_close(&bb, comm);
    }

    big_file_mpi_close(&bf, comm);

    /* return the values to the correct location */
    mpsort_mpi(buf, localsize, sizeof(buf[0]), _radix, 8, NULL, comm);

    /* read out the values */
    localsize = 0;
    for(pm_kiter_init(pm, &kiter);
       !pm_kiter_stop(&kiter);
        pm_kiter_next(&kiter)) {
        data[kiter.ind] = buf[localsize].value[0];
        data[kiter.ind + 1] = buf[localsize].value[1];
        localsize ++;
    }

    free(buf);
    return 0;
}
