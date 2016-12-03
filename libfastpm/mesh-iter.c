#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <fastpm/libfastpm.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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

