#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bigfile.h>
#include <bigfile-mpi.h>
#include <mpsort.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/string.h>
#include <fastpm/logging.h>

#include <fastpm/io.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

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

struct sort_data {
    FastPMStore * p; /* temporary store for decoding */
    FastPMPackingPlan * plan;
};

void
FastPMSnapshotSortByID(const void * ptr, void * radix, void * arg)
{
    struct sort_data * data = arg;

    FastPMStore * p = data->p;
    FastPMPackingPlan * plan = data->plan;

    fastpm_packing_plan_unpack_ci(plan, FASTPM_STORE_COLUMN_INDEX(id), p, 0, (void*) ptr);

    *((uint64_t*) radix) = p->id[0];
}

void
FastPMSnapshotSortByLength(const void * ptr, void * radix, void * arg)
{
    struct sort_data * data = arg;

    FastPMStore * p = data->p;
    FastPMPackingPlan * plan = data->plan;

    fastpm_packing_plan_unpack_ci(plan, FASTPM_STORE_COLUMN_INDEX(length), p, 0, (void*) ptr);

    *((uint64_t*) radix) = -p->length[0];
}

void
FastPMSnapshotSortByAEmit(const void * ptr, void * radix, void * arg)
{
    struct sort_data * data = arg;

    FastPMStore * p = data->p;
    FastPMPackingPlan * plan = data->plan;

    fastpm_packing_plan_unpack_ci(plan, FASTPM_STORE_COLUMN_INDEX(aemit), p, 0, (void*) ptr);

    /* larger than 53 is probably good enough but perhaps should have used ldexp (>GLIBC 2.19)*/
    *((uint64_t*) radix) = p->aemit[0] * (1L << 60L);
}

void
fastpm_sort_snapshot(FastPMStore * p, MPI_Comm comm, FastPMSnapshotSorter sorter, int redistribute)
{
    int64_t size = p->np;
    int NTask;
    int ThisTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);

    size_t localsize = size * (ThisTask + 1) / NTask - size * ThisTask / NTask;

    FastPMPackingPlan plan[1];
    fastpm_packing_plan_init(plan, p, p->attributes);
    size_t elsize = plan->elsize;

    if(redistribute) {
        if(MPIU_Any(comm, localsize > p->np_upper)) {
            fastpm_info("redistribution is requested in sorting, but the store is not large enough, will not redistribute: \n", p->np_upper, localsize);
            localsize = p->np;
        }
    } else {
        localsize = p->np;
    }

    void * send_buffer = malloc(elsize * p->np);
    void * recv_buffer = malloc(elsize * localsize);
    ptrdiff_t i;

    FastPMStore ptmp[1];
    fastpm_store_init(ptmp, 1, p->attributes, FASTPM_MEMORY_HEAP);

    for(i = 0; i < p->np; i ++) {
        fastpm_packing_plan_pack(plan, p, i, send_buffer + i * plan->elsize);
    }

    struct sort_data data[1];
    data->p = ptmp;
    data->plan = plan;

    mpsort_mpi_newarray(send_buffer, p->np, recv_buffer, localsize, elsize, sorter, 8, data, comm);

    for(i = 0; i < localsize; i ++) {
        fastpm_packing_plan_unpack(plan, p, i, recv_buffer + i * plan->elsize);
    }
    p->np = localsize;
    fastpm_store_destroy(ptmp);
    free(recv_buffer);
    free(send_buffer);

    double min, max, std, mean;

    MPIU_stats(comm, p->np, "<->s",
        &min, &mean, &max, &std);

    fastpm_info("number of particles after sorting: min = %g max = %g mean = %g std = %g\n",
            min, max, mean, std);
}

void
write_snapshot_header(FastPMSolver * fastpm, FastPMStore * p,
    BigFile * bf, MPI_Comm comm)
{
    BigBlock bb;
    if(0 != big_file_mpi_create_block(bf, &bb, "Header", "i8", 0, 1, 0, comm)) {
        fastpm_raise(-1, "Failed to create the header block: %s\n", big_file_get_error_message());
    }

    double H0 = 100.;
    /* Conversion from peculiar velocity to RSD,
     * http://mwhite.berkeley.edu/Talks/SantaFe12_RSD.pdf */
    double RSD = 1.0 / (H0 * p->a_x * HubbleEa(p->a_x, fastpm->cosmology));

    fastpm_info("RSD factor %e\n", RSD);

    double ScalingFactor = p->a_x;
    double OmegaM = fastpm->cosmology->OmegaM;
    double OmegaLambda = fastpm->cosmology->OmegaLambda;
    double HubbleParam = fastpm->config->hubble_param;
    double BoxSize = fastpm->config->boxsize;
    double rho_crit = 27.7455; /* 1e10 Msun /h*/
    uint64_t NC = fastpm->config->nc;

    double M0 = OmegaM * rho_crit * (BoxSize / NC) * (BoxSize / NC) * (BoxSize / NC);
    double MassTable[6] = {0, M0, 0, 0, 0, 0};
    uint64_t TotNumPart[6] = {0, fastpm_store_get_np_total(p, comm), 0, 0, 0, 0};

    /* FIXME: move some of these to fastpm.c; see if we can reduce the number of entries in fastpm->config. */
    big_block_set_attr(&bb, "NC", &NC, "i8", 1);
    big_block_set_attr(&bb, "BoxSize", &BoxSize, "f8", 1);
    big_block_set_attr(&bb, "ScalingFactor", &ScalingFactor, "f8", 1);
    big_block_set_attr(&bb, "RSDFactor", &RSD, "f8", 1);
    big_block_set_attr(&bb, "OmegaM", &OmegaM, "f8", 1);
    big_block_set_attr(&bb, "OmegaLambda", &OmegaLambda, "f8", 1);
    big_block_set_attr(&bb, "HubbleParam", &HubbleParam, "f8", 1);
    big_block_set_attr(&bb, "M0", &M0, "f8", 1);
    big_block_set_attr(&bb, "LibFastPMVersion", LIBFASTPM_VERSION, "S1", strlen(LIBFASTPM_VERSION));

    /* Compatibility with MP-Gadget */
    double UnitVelocity_in_cm_per_s = 1e5; /* 1 km/sec */
    double UnitLength_in_cm = 3.085678e21 * 1e3; /* 1.0 Mpc /h */
    double UnitMass_in_g = 1.989e43;       /* 1e10 Msun/h*/
    int UsePeculiarVelocity = 1;

    big_block_set_attr(&bb, "Omega0", &OmegaM, "f8", 1);
    big_block_set_attr(&bb, "TotNumPart", &TotNumPart, "i8", 6);
    big_block_set_attr(&bb, "MassTable", MassTable, "f8", 6);
    big_block_set_attr(&bb, "Time", &ScalingFactor, "f8", 1);
    big_block_set_attr(&bb, "UsePeculiarVelocity", &UsePeculiarVelocity, "i4", 1);
    big_block_set_attr(&bb, "UnitLength_in_cm", &UnitLength_in_cm, "f8", 1);
    big_block_set_attr(&bb, "UnitMass_in_g", &UnitMass_in_g, "f8", 1);
    big_block_set_attr(&bb, "UnitVelocity_in_cm_per_s", &UnitVelocity_in_cm_per_s, "f8", 1);
    big_block_mpi_close(&bb, comm);
}

int
write_snapshot_data(FastPMStore * p,
        int Nwriters,
        int append,
        BigFile * bf,
        const char * dataset,
        MPI_Comm comm
)
{
    int64_t size = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_LONG, MPI_SUM, comm);

    int NTask;
    MPI_Comm_size(comm, &NTask);

    /* at most 16 writers per file;
     * this would mean each writer writes about 2 M items and would trigger
     * aggregated IO when Nfile is very small. (items are typically less than 12 byes each.
     * */

    if(Nwriters == 0 || Nwriters > NTask) Nwriters = NTask;

    size_t items_per_file = 32 * 1024 * 1024;
    size_t writers_per_file = 16; /* at least use this many on small datasets */
    size_t max_elsize = 36; /* FIXME: find this from all columns */

    /* for smaller data sets, aggregate and write. (way faster) */
    /* 12 is a number to adjust for slight imbalance between nodes */
    big_file_mpi_set_aggregated_threshold((items_per_file / writers_per_file + 12) * max_elsize);

    int Nfile = (size + items_per_file - 1) / items_per_file;

    if(Nfile < 1) Nfile = 1;

    if(Nwriters > Nfile * writers_per_file) Nwriters = Nfile * writers_per_file;


    fastpm_info("Writing %ld objects to %d files with %d writers\n", size, Nfile, Nwriters);

    #define DEFINE_COLUMN_IO(name, dtype_, column) \
        {name, dtype_, \
                    FASTPM_STORE_COLUMN_INFO(p, column).dtype, \
                    FASTPM_STORE_COLUMN_INFO(p, column).nmemb, \
                    FASTPM_STORE_COLUMN_INFO(p, column).attribute, \
                    FASTPM_STORE_COLUMN_INDEX(column) \
                }

    struct {
        char * name;
        char * dtype_out;
        char * dtype;
        int nmemb;
        FastPMColumnTags attribute;
        int ci;
    } * descr, BLOCKS[] = {
        DEFINE_COLUMN_IO("Position",        "f4", x),
        DEFINE_COLUMN_IO("InitialPosition", "f4", q),
        DEFINE_COLUMN_IO("DX1",             "f4", dx1),
        DEFINE_COLUMN_IO("DX2",             "f4", dx2),
        DEFINE_COLUMN_IO("Velocity",        "f4", v),
        DEFINE_COLUMN_IO("ID",              "i8", id),
        DEFINE_COLUMN_IO("Aemit",           "f4", aemit),
        DEFINE_COLUMN_IO("Potential",       "f4", potential),
        DEFINE_COLUMN_IO("Density",         "f4", rho),
        DEFINE_COLUMN_IO("Tidal",           "f4", tidal),
        DEFINE_COLUMN_IO("Length",          "i4", length),
        DEFINE_COLUMN_IO("MinID",           "i8", minid),
        DEFINE_COLUMN_IO("Task",            "i4", task),
        DEFINE_COLUMN_IO("Rdisp",           "f4", rdisp),
        DEFINE_COLUMN_IO("Vdisp",           "f4", vdisp),
        DEFINE_COLUMN_IO("RVdisp",          "f4", rvdisp),
        {NULL, },
    };
    if (!append) {
        BigBlock bb;
        /* create the root block for the dataset attributes specific to this dataset. */
        if(0 != big_file_mpi_create_block(bf, &bb, dataset, NULL, 0, 0, 0, comm)) {
            fastpm_raise(-1, "Failed to create the block: %s\n", big_file_get_error_message());
        }
        big_block_mpi_close(&bb, comm);
    }
    for(descr = BLOCKS; descr->name; descr ++) {
        if(p->columns[descr->ci] == NULL) continue;

        fastpm_info("Writing block %s of (%s, %d)\n", descr->name, descr->dtype, descr->nmemb);
        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        char * blockname = fastpm_strdup_printf("%s/%s", dataset, descr->name);
        if(!append) {
            if(0 != big_file_mpi_create_block(bf, &bb, blockname, descr->dtype_out, descr->nmemb,
                        Nfile, size, comm)) {
                fastpm_raise(-1, "Failed to create the block: %s\n", big_file_get_error_message());
            }
            big_block_seek(&bb, &ptr, 0);
        } else {
            if(0 != big_file_mpi_open_block(bf, &bb, blockname, comm)) {
                /* if open failed, create an empty block instead.*/
                if(0 != big_file_mpi_create_block(bf, &bb, blockname, descr->dtype_out, descr->nmemb,
                            0, 0, comm)) {
                    fastpm_raise(-1, "Failed to create the block: %s\n", big_file_get_error_message());
                }
            }
            size_t oldsize = bb.size;
            /* FIXME : check the dtype and nmemb are consistent */
            big_file_mpi_grow_block(bf, &bb, Nfile, size, comm);
            big_block_seek(&bb, &ptr, oldsize);
        }

        /* packing the single column for IO */
        FastPMPackingPlan plan[1];
        fastpm_packing_plan_init(plan, p, descr->attribute);

        size_t elsize = plan->elsize;
        void * buffer = malloc(elsize * p->np);

        ptrdiff_t i;
        for(i = 0; i < p->np; i ++) {
            fastpm_packing_plan_pack(plan, p, i, ((char *) buffer) + i * elsize);
        }

        big_array_init(&array, buffer, descr->dtype, 2, (size_t[]) {p->np, descr->nmemb}, NULL );
        big_block_mpi_write(&bb, &ptr, &array, Nwriters, comm);
        free(buffer);

        big_block_mpi_close(&bb, comm);
        free(blockname);
    }

    return 0;
}

int
write_snapshot(FastPMSolver * fastpm, FastPMStore * p,
        const char * filebase,
        const char * dataset,
        int Nwriters)
{
    CLOCK(meta);

    int ThisTask = fastpm->ThisTask;
    MPI_Comm comm = fastpm->comm;

    ENTER(meta);
    if (ThisTask == 0)
        fastpm_path_ensure_dirname(filebase);

    LEAVE(meta);

    MPI_Barrier(comm);

    BigFile bf;
    if(0 != big_file_mpi_create(&bf, filebase, comm)) {
        fastpm_raise(-1, "Failed to create the file: %s\n", big_file_get_error_message());
    }

    write_snapshot_header(fastpm, p, &bf, comm);

    fastpm_info("Writring a catalog to %s [%s] with\n", filebase, dataset);

    write_snapshot_data(p,  Nwriters, 0, &bf, dataset, comm);

    big_file_mpi_close(&bf, comm);

    return 0;
}

int
append_snapshot(FastPMSolver * fastpm, FastPMStore * p,
        const char * filebase,
        const char * dataset,
        int Nwriters)
{
    CLOCK(meta);

    int ThisTask = fastpm->ThisTask;
    MPI_Comm comm = fastpm->comm;

    ENTER(meta);
    if (ThisTask == 0)
        fastpm_path_ensure_dirname(filebase);
    MPI_Barrier(comm);
    LEAVE(meta);


    BigFile bf;
    if(0 != big_file_mpi_open(&bf, filebase, comm)) {
        fastpm_raise(-1, "Failed to create the file: %s\n", big_file_get_error_message());
    }

    fastpm_info("Appending a catalog to %s [%s]\n", filebase, dataset);
    write_snapshot_data(p, Nwriters, 1, &bf, dataset, comm);

    big_file_mpi_close(&bf, comm);

    return 0;
}

int
read_snapshot(FastPMSolver * fastpm, FastPMStore * p, const char * filebase)
{
    return 0;
}

struct BufType {
    uint64_t ind;
    uint64_t ind2;
    float value[2];
};

static void
_radix(const void * ptr, void * radix, void * arg)
{
    const struct BufType * buf = ptr;
    uint64_t * key = radix;
    *key = buf->ind;
}

static void
_radix2(const void * ptr, void * radix, void * arg)
{
    const struct BufType * buf = ptr;
    uint64_t * key = radix;
    *key = buf->ind2;
}

int
write_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters)
{
    CLOCK(meta);

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

    ENTER(meta);

    if (ThisTask == 0)
        fastpm_path_ensure_dirname(filename);

    MPI_Barrier(comm);
    LEAVE(meta);

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

    if(0 != big_file_mpi_create(&bf, filename, comm)) {
        fastpm_raise(-1, "Failed to create the file: %s\n", big_file_get_error_message());
    }
    {
        BigBlock bb;
        BigArray array;
        BigBlockPtr ptr;
        if(0 != big_file_mpi_create_block(&bf, &bb, blockname, "c8", 1, Nfile, size, comm)) {
            fastpm_raise(-1, "Failed to create the block : %s\n", big_file_get_error_message());
        }
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
        if (0 != big_file_mpi_open_block(&bf, &bb, blockname, comm)) {
            fastpm_raise(-1, "Failed to open the block: %s\n", big_file_get_error_message());
        }

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

/**
 *
 * This routine creates an angular mesh from a data file that stores RA, DEC in degrees
 *
 * store : particle storage for the file; must be preallocated large enough to handle the file.
 * filename : location of bigfile data repository, must contain RA and DEC columns
 * r : an array of radial distance to place the mesh
 * aemit : an array of emission redshift; that is also stored in the mesh; This shall be solved from r
 * before the function is called.
 * Nr : number of entries in aemit / r
 * comm : MPI Communicator.
 *
 * The ->x and ->aemit elements are updated.
 * np is also updated.
 * XXX There may be issues at low-z as RA-DEC grid is concentrated around the observer. In that
 * case may need to downsample the RA-DEC grid points
 */
size_t
read_angular_grid(FastPMStore * store,
        const char * filename,
        const double * r,
        const double * aemit,
        const size_t Nr,
        int sampling_factor,
        MPI_Comm comm)
{

    int NTask;
    int ThisTask;

    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    BigFile bf;
    BigBlock bb;
    BigArray array;
    BigBlockPtr ptr;

    if(0 != big_file_mpi_open(&bf, filename, comm)) {
        goto exc_open;
    }

    if(0 != big_file_mpi_open_block(&bf, &bb, "RA", comm)) {
        goto exc_open;
    }
    size_t localstart = bb.size * ThisTask / NTask;
    size_t localend = bb.size * (ThisTask + 1) / NTask;
    size_t localsize = localend - localstart;

    fastpm_info("Reading %td ra dec coordinates with sampling of %td; localsize = %td\n",
                bb.size, sampling_factor, localsize);

    if(0 != big_block_mpi_close(&bb, comm)) {
        goto exc_open;
    }

    if (store==NULL){
      if(0 != big_file_mpi_close(&bf, comm)) {
          goto exc_open;
      }
      return localsize*Nr/sampling_factor+1;//XXX check for proper rounding off factor
    }


    if(0 != big_file_mpi_open_block(&bf, &bb, "RA", comm)) {
        goto exc_open_blk;
    }
    if(0 != big_block_seek(&bb, &ptr, 0)) {
        goto exc_seek;
    }

    double * RA = malloc(localsize * sizeof(double));
    double * DEC = malloc(localsize * sizeof(double));

    if(0 != big_array_init(&array, RA, "f8", 1, (size_t[]) {localsize, }, NULL)) {
        goto exc_arr;
    }
    if(0 != big_block_mpi_read(&bb, &ptr, &array, NTask, comm)) {
        goto exc_read;
    }
    if(0 != big_block_mpi_close(&bb, comm)) {
        goto exc_close_blk;
    }

    if(0 != big_file_mpi_open_block(&bf, &bb, "DEC", comm)) {
        goto exc_open_blk;
    }
    if(0 != big_block_seek(&bb, &ptr, 0)) {
        goto exc_seek;
    }
    if(0 != big_array_init(&array, DEC, "f8", 1, (size_t[]) {localsize, }, NULL)) {
        goto exc_arr;
    }
    if(0 != big_block_mpi_read(&bb, &ptr, &array, NTask, comm)) {
        goto exc_read;
    }
    if(0 != big_block_mpi_close(&bb, comm)) {
        goto exc_close_blk;
    }

    if(0 != big_file_mpi_close(&bf, comm)) {
        goto exc_close;
    }


    double * x = malloc(localsize * sizeof(double));
    double * y = malloc(localsize * sizeof(double));
    double * z = malloc(localsize * sizeof(double));



    ptrdiff_t i;
    ptrdiff_t j;
    ptrdiff_t n;
    n = store->np;
    double d2r=180./M_PI;
    for(i = 0; i < localsize; i ++) {
        RA[i] /= d2r;
        //DEC[i] /= d2r;
        DEC[i]=M_PI/2.-DEC[i]/d2r;
        /* FIXME conversion is likely wrong. */
        x[i] = sin(DEC[i]) * cos(RA[i]);
        y[i] = sin(DEC[i]) * sin(RA[i]);
        z[i] = cos(DEC[i]);
    }

    for(j = 0; j < Nr; j ++) {
      for(i = 0; i < localsize; i+=sampling_factor) {
            store->x[n][0] = x[i] * r[j];
            store->x[n][1] = y[i] * r[j];
            store->x[n][2] = z[i] * r[j];
            store->aemit[n] = aemit[j];
            n++;
        }
        if(n == store->np_upper) {
            fastpm_raise(-1, "Too many grid points on the grid, the limit is %td with %td r bins, i=%td  j=%td n=%td \n", store->np_upper,Nr,i,j,n);
        }
    }

    store->np = n;

    fastpm_info("Generated %td x, y, z coordinates locally\n", n);

    free(DEC);
    free(RA);
    free(x);
    free(y);
    free(z);

    return n;

    exc_close:
    exc_close_blk:
    exc_read:
    exc_arr:
        free(DEC);
        free(RA);
    exc_seek:
    exc_open_blk:
    exc_open:

        fastpm_raise(-1, "Failed to read angular file %s, for %s\n", filename, big_file_get_error_message());

    return 0;
}

void
write_snapshot_attr(const char * filebase,
    const char * dataset,
    const char * attrname,
    void * buf,
    const char * dtype,
    size_t nmemb,
    MPI_Comm comm)
{
    BigFile bf;
    BigBlock bb;
    if(0 != big_file_mpi_open(&bf, filebase, comm)) {
        fastpm_raise(-1, "Failed to open the file: %s\n", big_file_get_error_message());
    }

    if(0 != big_file_mpi_open_block(&bf, &bb, dataset, comm)) {
        fastpm_raise(-1, "Failed to open the dataset : %s\n", big_file_get_error_message());
    }

    big_block_set_attr(&bb, attrname, buf, dtype, nmemb);

    big_block_mpi_close(&bb, comm);
    big_file_mpi_close(&bf, comm);
}

void
write_aemit_hist(const char * filebase, const char * dataset,
            int64_t * hist,
            double * aedges,
            size_t nedges,
            MPI_Comm comm)
{
    /* write an index to the dataset (usually Header block) */
    /* offset[i]: offset[i] + size[i] contains all particles from aemit[i]: aemit[i+1] */
    BigFile bf;
    BigBlock bb;

    if(0 != big_file_mpi_open(&bf, filebase, comm)) {
        fastpm_raise(-1, "Failed to open the file: %s\n", big_file_get_error_message());
    }

    if(0 != big_file_mpi_open_block(&bf, &bb, dataset, comm)) {
        fastpm_raise(-1, "Failed to open the dataset : %s\n", big_file_get_error_message());
    }

    /* starting edges of the bins */
    big_block_remove_attr(&bb, "aemitIndex.edges");
    big_block_set_attr(&bb, "aemitIndex.edges", aedges, "f8", nedges);

    /* number in each layer bin 0 is the first layer, since it is the only layer outside of the edges */
    big_block_remove_attr(&bb, "aemitIndex.size");
    big_block_set_attr(&bb, "aemitIndex.size", hist, "i8", nedges + 1);

    int64_t * offset = malloc(sizeof(int64_t) * (nedges + 1));
    int64_t * size = malloc(sizeof(int64_t) * nedges);

    /* the starting particles offset for each layer */
    offset[0] = 0;
    int i;
    for(i = 1; i < nedges + 1; i ++) {
        offset[i] = offset[i - 1] + hist[i];
    }
    big_block_remove_attr(&bb, "aemitIndex.offset");
    big_block_set_attr(&bb, "aemitIndex.offset", offset, "i8", nedges + 1);

    big_block_mpi_close(&bb, comm);
    big_file_mpi_close(&bf, comm);

    free(offset);
    free(size);
}

