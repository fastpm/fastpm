FASTPM_BEGIN_DECLS

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
