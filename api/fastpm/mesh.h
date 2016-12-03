
FASTPM_BEGIN_DECLS

typedef struct {
    /* in units of real numbers, not bytes. */
    int ndim;
    ptrdiff_t start[3];
    ptrdiff_t shape[3];
    ptrdiff_t strides[3];
    ptrdiff_t size;
} FastPMArrayLayout;

typedef struct FastPMMeshPrivate FastPMMeshPrivate;
typedef struct {
    FastPMMemory *mem;
    int ndim;

    FastPMArrayLayout ral;
    FastPMArrayLayout cal;

    int Nproc[2];
    int NTask;
    double BoxSize[3];
    ptrdiff_t Nmesh[3];

    double CellSize[3];
    double InvCellSize[3];
    double Volume;
    double Norm;
    FastPMMeshPrivate * priv;
    MPI_Comm comm;
} FastPMMesh;


void
fastpm_mesh_init(FastPMMesh * self,
        int ndim,
        double BoxSize,
        ptrdiff_t Nmesh,
        ptrdiff_t Nproc[],
        MPI_Comm comm);

void
fastpm_mesh_destroy(FastPMMesh * self);

FastPMFloat *
fastpm_mesh_alloc(FastPMMesh * self);

void fastpm_mesh_free(FastPMMesh * self, FastPMFloat * array);

void
fastpm_mesh_r2c(FastPMMesh * self, FastPMFloat * r, FastPMFloat * c);

void
fastpm_mesh_c2r(FastPMMesh * self, FastPMFloat * c, FastPMFloat * r);

void
fastpm_mesh_copy(FastPMMesh * self, FastPMFloat * from, FastPMFloat * to);

int
fastpm_mesh_ipos_to_rank(FastPMMesh * self, int ipos[]);

int
fastpm_mesh_pos_to_ipos(FastPMMesh * self, double pos, int dim);

FASTPM_END_DECLS
