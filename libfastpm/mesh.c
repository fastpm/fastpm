#include <string.h>
#include <alloca.h>
#include <stdlib.h>

#include <mpi.h>
#include <pfft.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#if FASTPM_FFT_PRECISION == 64
#define _pfft_plan pfft_plan
#define _pfft_create_procmesh pfft_create_procmesh
#define _pfft_local_size_dft_r2c pfft_local_size_dft_r2c
#define _pfft_plan_dft_r2c pfft_plan_dft_r2c
#define _pfft_plan_dft_c2r pfft_plan_dft_c2r
#define _pfft_execute_dft_r2c pfft_execute_dft_r2c
#define _pfft_execute_dft_c2r pfft_execute_dft_c2r
#define _pfft_init pfft_init
#define _pfft_cleanup pfft_cleanup
#define _pfft_destroy_plan pfft_destroy_plan
#else
#define _pfft_plan pfftf_plan
#define _pfft_create_procmesh pfftf_create_procmesh
#define _pfft_local_size_dft_r2c pfftf_local_size_dft_r2c
#define _pfft_plan_dft_r2c pfftf_plan_dft_r2c
#define _pfft_plan_dft_c2r pfftf_plan_dft_c2r
#define _pfft_execute_dft_r2c pfftf_execute_dft_r2c
#define _pfft_execute_dft_c2r pfftf_execute_dft_c2r
#define _pfft_init pfftf_init
#define _pfft_cleanup pfftf_cleanup
#define _pfft_destroy_plan pfftf_destroy_plan
#endif

struct FastPMMeshPrivate {
    MPI_Datatype ptrtype;
    MPI_Comm cart;
    int * Mesh2Rank[2];
    size_t allocsize;
    _pfft_plan plan_r2c;
    _pfft_plan plan_r2c_inplace;
    _pfft_plan plan_c2r;
    _pfft_plan plan_c2r_inplace;
};

void
fastpm_mesh_init(FastPMMesh * self,
        double BoxSize,
        ptrdiff_t Nmesh,
        int ndim,
        ptrdiff_t Nproc[],
        MPI_Comm comm)
{
    self->mem = _libfastpm_get_gmem();
    self->priv = malloc(sizeof(self->priv[0]));
    FastPMMeshPrivate * priv = self->priv;
    if(sizeof(ptrdiff_t) == 8) {
        priv->ptrtype = MPI_LONG;
    } else {
        priv->ptrtype = MPI_INT;
    }

    self->ndim = ndim;
    self->ral.ndim = ndim;
    self->cal.ndim = ndim;
    int d = 0;
    self->Volume = 1;
    self->Norm = 1;
    for(d = 0; d < ndim; d ++) {
        self->BoxSize[d] = BoxSize;
        self->Nmesh[d] = Nmesh;
        self->CellSize[d] = BoxSize / Nmesh;
        self->InvCellSize[d] = Nmesh / BoxSize;
        self->Volume *= BoxSize;
        self->Norm *= Nmesh;
    }
    self->NTask = 1;
    for(d = 0; d < ndim - 1; d ++) {
        self->Nproc[d] = Nproc[d];
        self->NTask *= Nproc[d];
    }
    self->comm = comm;
    int NTask;
    MPI_Comm_size(comm, &NTask);

    if(NTask != self->NTask) {
        fastpm_raise(-1, "NTask and Nproc mismatch\n");
    }

    _pfft_create_procmesh(2, comm, self->Nproc, &priv->cart);

    priv->allocsize = 2 * _pfft_local_size_dft_r2c(
            self->ndim, self->Nmesh, priv->cart, 
            PFFT_TRANSPOSED_OUT | PFFT_PADDED_R2C, 
            self->ral.shape, self->ral.start,
            self->cal.shape, self->cal.start);

    if(self->ndim != 3) {
        fastpm_raise(-1, "only three dimensions are supported for now\n");
    } else {
        self->ral.strides[2] = 1;
        self->ral.strides[1] = self->ral.shape[2];
        self->ral.strides[0] = self->ral.shape[1] * self->ral.strides[1];
        self->ral.size = self->ral.shape[0] * self->ral.strides[0];

        /* remove padding from the view */
        self->ral.shape[2] = self->Nmesh[2];

        /* PFFT transposed, y, z, x */
        self->cal.strides[0] = 1;
        self->cal.strides[2] = self->cal.shape[0];
        self->cal.strides[1] = self->cal.shape[2] * self->cal.strides[2];
        self->cal.size = self->cal.shape[1] * self->cal.strides[1];
    }

    for(d = 0; d < self->ndim - 1; d ++) {
        MPI_Comm projected;
        int remain_dims[100];
        int j;
        for(j = 0; j < self->ndim - 1; j ++) {
            remain_dims[j] = 0;
        }
        remain_dims[d] = 1; 
        priv->Mesh2Rank[d] = malloc(sizeof(ptrdiff_t) * self->Nmesh[d]);
        int * edges_int = (int*) malloc(sizeof(edges_int[0]) * (self->Nproc[d] + 1));

        MPI_Cart_sub(priv->cart, remain_dims, &projected);
        MPI_Allgather(&self->ral.start[d], 1, priv->ptrtype, 
            edges_int, 1, priv->ptrtype, projected);
        int ntask;
        MPI_Comm_size(projected, &ntask);

        MPI_Comm_free(&projected);
        /* Last edge is at the edge of the box */
        edges_int[j] = self->Nmesh[d];
        /* fill in the look up table */
        for(j = 0; j < self->Nproc[d]; j ++) {
            int i;
            for(i = edges_int[j]; i < edges_int[j+1]; i ++) {
                priv->Mesh2Rank[d][i] = j;
            }
        }
        free(edges_int);
    }

    FastPMFloat * workspace = fastpm_mesh_alloc(self);
    FastPMFloat * canvas = fastpm_mesh_alloc(self);
    priv->plan_r2c = _pfft_plan_dft_r2c(
            self->ndim, self->Nmesh, (void*) workspace, (void*) canvas, 
            priv->cart,
            PFFT_FORWARD, 
            PFFT_TRANSPOSED_OUT
            | PFFT_PADDED_R2C 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_r2c_inplace = _pfft_plan_dft_r2c(
            self->ndim, self->Nmesh, (void*) workspace, (void*) workspace, 
            priv->cart,
            PFFT_FORWARD, 
            PFFT_TRANSPOSED_OUT
            | PFFT_PADDED_R2C 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_c2r = _pfft_plan_dft_c2r(
            self->ndim, self->Nmesh, (void*) workspace, (void*) canvas, 
            priv->cart,
            PFFT_BACKWARD, 
            PFFT_TRANSPOSED_IN
            | PFFT_PADDED_C2R 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_c2r_inplace = _pfft_plan_dft_c2r(
            self->ndim, self->Nmesh, (void*) workspace, (void*) workspace, 
            priv->cart,
            PFFT_BACKWARD, 
            PFFT_TRANSPOSED_IN
            | PFFT_PADDED_C2R 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    fastpm_mesh_free(self, canvas);
    fastpm_mesh_free(self, workspace);
}

void
fastpm_mesh_destroy(FastPMMesh * self)
{
    int d;
    for(d = 0; d < self->ndim - 1; d ++) {
        free(self->priv->Mesh2Rank[d]);
    }
    MPI_Comm_free(&self->priv->cart);
    _pfft_destroy_plan(self->priv->plan_r2c);
    _pfft_destroy_plan(self->priv->plan_c2r);
    _pfft_destroy_plan(self->priv->plan_r2c_inplace);
    _pfft_destroy_plan(self->priv->plan_c2r_inplace);
    free(self->priv);
}

FastPMFloat *
fastpm_mesh_alloc(FastPMMesh * self)
{
    void * p = fastpm_memory_alloc(self->mem, sizeof(FastPMFloat) * self->priv->allocsize, FASTPM_MEMORY_HEAP);
    fastpm_memory_tag(self->mem, p, "FastPMFloat");
    return p;
}

void
fastpm_mesh_free(FastPMMesh * self, FastPMFloat * array)
{
    fastpm_memory_free(self->mem, array);
}

void
fastpm_mesh_r2c(FastPMMesh * self, FastPMFloat * r, FastPMFloat * c)
{
    if(c == r) {
        pfft_execute_dft_r2c(self->priv->plan_r2c_inplace, (void*)r, (void*)c);
    } else {
        pfft_execute_dft_r2c(self->priv->plan_r2c, (void*) r, (void*) c);
    }
}

void
fastpm_mesh_c2r(FastPMMesh * self, FastPMFloat * c, FastPMFloat * r)
{
    if(c == r) {
        pfft_execute_dft_c2r(self->priv->plan_c2r_inplace, (void*) r, (void*) c);
    } else {
        pfft_execute_dft_c2r(self->priv->plan_c2r, (void*) r, (void*) c);
    }
}

void
fastpm_mesh_copy(FastPMMesh * self, FastPMFloat * from, FastPMFloat * to)
{
    if(from == to) {
        return;
    } else {
        memcpy(to, from, sizeof(FastPMFloat) * self->priv->allocsize);
    }
}

int
fastpm_mesh_ipos_to_rank(FastPMMesh * self, int ipos, int dim)
{
    while(ipos >= self->Nmesh[dim]) ipos -= self->Nmesh[dim];
    while(ipos < 0) ipos += self->Nmesh[dim];

    return self->priv->Mesh2Rank[dim][ipos];
}

int
fastpm_mesh_pos_to_ipos(FastPMMesh * self, double pos, int dim)
{
    return (int) floor(self->InvCellSize[dim] * pos);
}


static int
_fastpm_mesh_iter_next_k(FastPMMeshIter * iter);
static int
_fastpm_mesh_iter_next_x(FastPMMeshIter * iter);
static float **
_fastpm_mesh_iter_make_table_k(
    FastPMMeshIter * iter,
    double (*func)(int i, double dx, int nmesh));


FastPMMeshIter *
fastpm_mesh_iter(FastPMMesh * self, enum FastPMMeshIterType type)
{
    FastPMMeshIter * iter = (FastPMMeshIter*) malloc(sizeof(iter[0]));
    iter->type = type;
#ifdef _OPENMP
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();
#else
    int nth = 1;
    int ith = 0;
#endif
    FastPMArrayLayout * al;
    switch(type) {
        case FASTPM_MESH_ITER_K:
            al = &self->cal;
        break;
        case FASTPM_MESH_ITER_X:
            al = &self->ral;
        break;
    }
    iter->start = ith * al->size / nth;
    iter->end = (ith + 1) * al->size / nth;
    int d;
    for(d = 0; d < self->ndim; d ++) {
        if(al->shape[d] == 0)
            iter->i[d] = 0;
        else
            iter->i[d] = (iter->start / al->strides[d]) % al->shape[d];

        iter->iabs[d] = iter->i[d] + al->start[d];
    }

    iter->mesh = self;

    switch(type) {
        case FASTPM_MESH_ITER_K:
            iter->start *= 2;
            iter->end *= 2;
            iter->next = _fastpm_mesh_iter_next_k;
            iter->table = _fastpm_mesh_iter_make_table_k;
        break;
        case FASTPM_MESH_ITER_X:
            iter->next = _fastpm_mesh_iter_next_x;
            iter->table = NULL;
        break;
    }
    iter->ind = iter->start;
    iter->al = al;
    iter->started = 0;
    return iter;
}

void
fastpm_mesh_iter_free(FastPMMeshIter * iter)
{
    free(iter);
}

int
fastpm_mesh_iter_next(FastPMMeshIter * iter)
{
    if(!iter->started) {
        iter->started = 1;
        return 1;
    }
    return iter->next(iter);
}

float **
fastpm_mesh_iter_make_table(
    FastPMMeshIter * iter,
    double (*func)(int ii, double dx, int nmesh))
{
    return iter->table(iter, func);
}

static int
_fastpm_mesh_iter_next_x(FastPMMeshIter * iter)
{
    iter->ind ++;
    int d;
    int carry = 1;
    for(d = iter->al->ndim - 1 ; carry && d >= 0; d--) {
        iter->i[d] += carry;
        carry = 0;
        if(iter->i[d] == iter->al->shape[d]) {
            iter->i[d] = 0;
            if(d == iter->al->ndim - 1) {
                /* Special padding for real array, ind shall increase two more */
                iter->ind += 2;
            }
            carry = 1;
        }
    }
    for(d = 0; d < iter->al->ndim; d ++)
        iter->iabs[d] = iter->i[d] + iter->al->start[d];
    return carry == 0;
}

static float **
_fastpm_mesh_iter_make_table_k(
    FastPMMeshIter * iter,
    double (*func)(int i, double dx, int nmesh))
{
    int n = 0;
    int d;
    for(d = 0; d < iter->mesh->ndim; d++) {
        n += iter->mesh->Nmesh[d];
    }
    char * buf = malloc(sizeof(float*) * iter->mesh->ndim + sizeof(float) * n);

    float ** rt = (float**) buf;
    buf += sizeof(float*) * iter->mesh->ndim;
    rt[0] = (float*) (buf);

    for(d = 1; d < iter->mesh->ndim; d++) {
        rt[d] = rt[d - 1] + iter->mesh->Nmesh[d];
    }

    for(d = 0; d < iter->mesh->ndim; d++) {
        int i;
        for(i = 0; i < iter->mesh->Nmesh[d]; i ++) {
            int ii = i;
            if(ii > iter->mesh->Nmesh[d] / 2) ii -= iter->mesh->Nmesh[d];
            rt[d][i] = func(ii, iter->mesh->CellSize[d], iter->mesh->Nmesh[d]);
        }
    }
    return rt;
}

static int
_fastpm_mesh_iter_next_k(FastPMMeshIter * iter)
{
    iter->ind += 2;
    /* transposed, y, z, x */
    int dd[] = {1, 2, 0};
    int d;
    int carry = 1;
    for(d = iter->al->ndim - 1 ; carry && d >= 0; d--) {
        iter->i[dd[d]] += carry;
        carry = 0;
        if(iter->i[dd[d]] == iter->al->shape[dd[d]]) {
            iter->i[dd[d]] = 0;
            carry = 1;
        }
    }
    for(d = 0; d < iter->al->ndim; d ++)
        iter->iabs[d] = iter->i[d] + iter->al->start[d];
    return carry == 0;
}

double
FastPMMeshK(int ii, double dx, int nmesh)
{
    if(ii >= nmesh / 2) return 0;
    double k = 2 * M_PI / (dx * nmesh) * ii;
    return k;
}

double
FastPMMeshKK(int ii, double dx, int nmesh)
{
    double k = 2 * M_PI / (dx * nmesh) * ii;

    return k * k;
}

static double
sinxoverx(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}

double
FastPMMeshKK3Point(int ii, double dx, int nmesh)
{
    double k = 2 * M_PI / (dx * nmesh) * ii;
    double w = k * dx;

    float ff1 = sinxoverx(0.5 * w);

    return k * k * ff1 * ff1;
}

double
FastPMMeshKK5Point(int ii, double dx, int nmesh)
{
    double k = 2 * M_PI / (dx * nmesh) * ii;
    double w = k * dx;

    float ff1 = sinxoverx(0.5 * w);
    float ff2 = sinxoverx(w);

    return k * k * (4 / 3.0 * ff1 * ff1 - 1 / 3.0 * ff2 * ff2);
}

double
FastPMMeshK4Point(int ii, double dx, int nmesh)
{

    if(ii >= nmesh / 2) return 0;
    double w;

    w = 2 * M_PI / nmesh * ii;

    /* order N = 1 super lanzcos kernel */
    /* 
     * This is the same as GADGET-2 but in fourier space: 
     * see gadget-2 paper and Hamming's book.
     * c1 = 2 / 3, c2 = 1 / 12
     * */
    return 1 / dx * 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}

