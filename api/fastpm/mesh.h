
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
fastpm_mesh_ipos_to_rank(FastPMMesh * self, int ipos, int dim);

int
fastpm_mesh_pos_to_ipos(FastPMMesh * self, double pos, int dim);

enum FastPMMeshIterType
{FASTPM_MESH_ITER_X, FASTPM_MESH_ITER_K };

typedef struct FastPMMeshIter FastPMMeshIter;

struct FastPMMeshIter {
    enum FastPMMeshIterType type;
    FastPMArrayLayout * al;
    ptrdiff_t start;
    ptrdiff_t end;
    ptrdiff_t ind;

    ptrdiff_t i[32];
    ptrdiff_t iabs[32];

    FastPMMesh * mesh;
    int (*next) (FastPMMeshIter * iter);
    float ** (*table) (FastPMMeshIter * iter, double (*func)(int ii, double dx, int nmesh));
    int started;
};

FastPMMeshIter *
fastpm_mesh_iter(FastPMMesh * self, enum FastPMMeshIterType type);

void
fastpm_mesh_iter_free(FastPMMeshIter * iter);

int
fastpm_mesh_iter_next(FastPMMeshIter * iter);

float **
fastpm_mesh_iter_make_table(
    FastPMMeshIter * iter,
    double (*func)(int ii, double dx, int nmesh));

double
FastPMMeshK(int ii, double dx, int nmesh);

double
FastPMMeshK4Point(int ii, double dx, int nmesh);

double
FastPMMeshKK3Point(int ii, double dx, int nmesh);

double
FastPMMeshKK(int ii, double dx, int nmesh);

double
FastPMMeshKK5Point(int ii, double dx, int nmesh);

FASTPM_END_DECLS
