#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"
#include "pmghosts.h"

/* from cic.c */
void fastpm_painter_init_cic(FastPMPainter * painter);

static void
_generic_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight, int diffdir);
static double
_generic_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], int diffdir);

static double
_linear_kernel(double x, double invh) {
    return 1.0 - fabs(x * invh);
}

static double
_linear_diff(double x, double invh) {
    if( x < 0) {
        return 1 * invh;
    } else {
        return - 1 * invh;
    }
}

static double
_quad_kernel(double x, double invh) {
    /*
     * Take from https://arxiv.org/abs/0804.0070
     * */
    x = fabs(x) * invh;
    if(x <= 0.5) {
        return 0.75 - x * x;
    } else {
        x = 1.5 - x;
        return (x * x) * 0.5;
    }
}

static double
_quad_diff(double x, double invh) {
    double factor;
    x *= invh;
    if ( x < 0) {
        x = -x;
        factor = -1 * invh;
    } else {
        factor = +1 * invh;
    }

    if(x < 0.5) {
        return factor * (- 2 * x);
    } else {
        return factor * (- (1.5 - x));
    }
}


static inline double __cached__(int *status, double * table, double x, double (*func)(double)){
    const double dx = 1e-3;
    const double tablemax = dx * 16384;
    const double tablemin = dx * 1;
    if(!*status) {
        int i;
        for(i = 0; i < 16384; i ++) {
            double x = dx * i;
            table[i] = func(x);
        }
        *status = 1;
    }
    if(x > tablemin && x < tablemax) {
        int i = fabs(x) / dx;
        return table[i];
    }
    return func(x);
}


static double __sinc__(double x) {
    x *= 3.1415927;
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}

static double __dsinc__(double x) {
    x *= 3.1415927;
    double r = 3.1415927;
    if(x < 1e-5 && x > -1e-5) {
        double xx = x * x;
        double xxxx = xx * xx;
        r *= - x / 3 + x*xx / 30 - xxxx*x/ 840 + xxxx * xx * x / 45360;
    } else {
        r *= 1 / x * cos(x) - 1 / (x *x) * sin(x);
    }
    return r;
}

static double
_lanczos_kernel(double x, double invh) {
    static int status = 0;
    static double table[16384];
    double s1 = __cached__(&status, table, x, __sinc__);
    double s2 = __cached__(&status, table, x * invh, __sinc__);
    return s1 * s2;
}

static double
_lanczos_diff(double x, double invh) {
    static int status = 0;
    static double table[16384];
    double u1 = __cached__(&status, table, x, __sinc__);
    double u2 = __cached__(&status, table, x, __dsinc__);
    double v1 = __cached__(&status, table, x * invh, __sinc__);
    double v2 = __cached__(&status, table, x * invh, __dsinc__) * invh;
    return u1 * v2 + u2 * v1;
}

void
fastpm_painter_init(FastPMPainter * painter, PM * pm,
    FastPMPainterType type, int support)
{
    painter->pm = pm;
    painter->paint = _generic_paint;
    painter->readout = _generic_readout;

    switch(type) {
        case FASTPM_PAINTER_CIC:
            fastpm_painter_init_cic(painter);
            painter->kernel = NULL;
            painter->diff = NULL;
        break;
        case FASTPM_PAINTER_LINEAR:
            painter->kernel = _linear_kernel;
            painter->diff = _linear_diff;
            support = 2;
        break;
        case FASTPM_PAINTER_QUAD:
            painter->kernel = _quad_kernel;
            painter->diff = _quad_diff;
            support = 3;
        break;
        case FASTPM_PAINTER_LANCZOS:
            painter->kernel = _lanczos_kernel;
            painter->diff = _lanczos_diff;
        break;
    }

    painter->support = support;
    painter->hsupport = 0.5 * support;
    painter->invh= 1 / (0.5 * support);
    painter->left = (support  - 1) / 2;
    painter->diffdir = -1;
    int nmax = 1;
    int d;
    for(d = 0; d < 3; d++) {
        nmax *= (support);
    }
    painter->Npoints = nmax;

    if (painter->support % 2 == 0){
        painter->shift = 0;
    } else {
        painter->shift = 0.5;
    }
}

void
fastpm_painter_init_diff(FastPMPainter * painter, FastPMPainter * base, int diffdir)
{
    *painter = *base;
    painter->diffdir = diffdir;
}

void
fastpm_painter_destroy(FastPMPainter * painter)
{
    /* nothing to do. */
    return;
}
static void
_fill_k(FastPMPainter * painter, double pos[3], int ipos[3], double k[3][64], int diffdir)
{
    PM * pm = painter->pm;
    double gpos[3];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->InvCellSize[d];
        ipos[d] = floor(gpos[d] + painter->shift) - painter->left;
        double dx = gpos[d] - ipos[d];
        int i;
        double sum = 0;
        for(i = 0; i < painter->support; i ++) {
            k[d][i] = painter->kernel(dx - i, painter->invh);
            sum += k[d][i];

            /*
             * norm is still from the true kernel,
             * but we replace the value with the derivative
             * */
            if(d == diffdir) {
                k[d][i] = painter->diff(dx - i, painter->invh) * pm->InvCellSize[d];
            }
        }
        /* normalize the kernel to conserve mass */
        for(i = 0; i < painter->support; i ++) {
                k[d][i] /= sum;
        }
        ipos[d] -= pm->IRegion.start[d];
    }
}

static void
_generic_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight, int diffdir)
{
    PM * pm = painter->pm;
    int ipos[3];
    /* the max support is 32 */
    double k[3][64];

    _fill_k(painter, pos, ipos, k, diffdir);

    int rel[3] = {0, 0, 0};
    int s2 = painter->support;
    while(rel[0] != s2) {
        double kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int r = rel[d];
            int targetpos = ipos[d] + r;
            kernel *= k[d][r];
            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(UNLIKELY(targetpos >= pm->IRegion.size[d]))
                goto outside;
            if(UNLIKELY(targetpos < 0))
                goto outside;
            ind += pm->IRegion.strides[d] * targetpos;
        }
#pragma omp atomic
        canvas[ind] += weight * kernel;

    outside:
        rel[2] ++;
        if(UNLIKELY(rel[2] == s2)) {
            rel[1] ++;
            rel[2] = 0;
        }
        if(UNLIKELY(rel[1] == s2)) {
            rel[1] = 0;
            rel[0] ++;
        }
        continue;
    }
    return;
}

static double
_generic_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], int diffdir)
{
    PM * pm = painter->pm;
    double value = 0;
    int ipos[3];
    double k[3][64];

    _fill_k(painter, pos, ipos, k, diffdir);

    int rel[3] = {0, 0, 0};

    int s2 = painter->support;
    while(rel[0] != s2) {
        double kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int r = rel[d];

            kernel *= k[d][r];

            int targetpos = ipos[d] + r;

            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(UNLIKELY(targetpos >= pm->IRegion.size[d])) {
                goto outside;
            }
            if(UNLIKELY(targetpos < 0))
                goto outside;
            ind += pm->IRegion.strides[d] * targetpos;
        }
        value += kernel * canvas[ind];
outside:
        rel[2] ++;
        if(UNLIKELY(rel[2] == s2)) {
            rel[1] ++;
            rel[2] = 0;
        }
        if(UNLIKELY(rel[1] == s2)) {
            rel[1] = 0;
            rel[0] ++;
        }
        continue;
    }
    return value;
}

void
fastpm_paint_local(FastPMPainter * painter, FastPMFloat * canvas,
    FastPMStore * p, size_t size,
    fastpm_posfunc get_position, int attribute)
{
    ptrdiff_t i;

    memset(canvas, 0, sizeof(canvas[0]) * painter->pm->allocsize);

    if(get_position == NULL) {
        get_position = p->get_position;
    }

#pragma omp parallel for
    for (i = 0; i < size; i ++) {
        double pos[3];
        double weight = attribute? p->to_double(p, i, attribute): 1.0;
        get_position(p, i, pos);
        painter->paint(painter, canvas, pos, weight, painter->diffdir);
    }
}

void
fastpm_paint(FastPMPainter * painter, FastPMFloat * canvas,
    FastPMStore * p, fastpm_posfunc get_position, int attribute)
{
    if(get_position == NULL) {
        get_position = p->get_position;
    }

    PMGhostData * pgd = pm_ghosts_create(painter->pm, p, p->attributes, get_position);

    fastpm_paint_local(painter, canvas, p, p->np + pgd->nghosts, get_position, attribute);

    pm_ghosts_free(pgd);
}

void
fastpm_readout_local(FastPMPainter * painter, FastPMFloat * canvas,
    FastPMStore * p, size_t size,
    fastpm_posfunc get_position, int attribute)
{

    ptrdiff_t i;
    if(get_position == NULL) {
        get_position = p->get_position;
    }
#pragma omp parallel for
    for (i = 0; i < size; i ++) {
        double pos[3];
        get_position(p, i, pos);
        double weight = painter->readout(painter, canvas, pos, painter->diffdir);
        p->from_double(p, i, attribute, weight);
    }
}

void
fastpm_readout(FastPMPainter * painter, FastPMFloat * canvas,
    FastPMStore * p, fastpm_posfunc get_position, int attribute)
{
    if(get_position == NULL) {
        get_position = p->get_position;
    }

    PMGhostData * pgd = pm_ghosts_create(painter->pm, p, p->attributes, get_position);

    fastpm_readout_local(painter, canvas, p, p->np + pgd->nghosts, get_position, attribute);

    pm_ghosts_reduce(pgd, attribute);
    pm_ghosts_free(pgd);
}
typedef struct {
    /* get_position is from base */
    /* from_double and to_double uses the pointers */
    FastPMStore base;
    FastPMStore * write_to;
    FastPMStore * read_from;
} DualStore;

static double
dual_store_to_double(FastPMStore * store, ptrdiff_t i, int attribute)
{
    DualStore * self = (DualStore*) store;

    return self->read_from->to_double(self->read_from, i, attribute);
}

static void
dual_store_from_double(FastPMStore * store, ptrdiff_t i, int attribute, double value)
{
    DualStore * self = (DualStore*) store;

    self->write_to->from_double(self->write_to, i, attribute, value);
}

static void
dual_store_init(DualStore * ds, FastPMStore * base)
{
    ds->base = * base;
    ds->base.to_double = dual_store_to_double;
    ds->base.from_double = dual_store_from_double;
}

static void
dual_store_set_write_to(DualStore * ds, FastPMStore * p)
{
    ds->write_to = p;
}

static void
dual_store_set_read_from(DualStore * ds, FastPMStore * p)
{
    ds->read_from = p;
}

void
fastpm_paint_gradient(FastPMPainter * painter, FastPMFloat * y,
    FastPMStore * p, fastpm_posfunc get_position, int attribute, FastPMStore * grad_pos, FastPMStore * grad_attr)
{
    DualStore dualstore[1];

    dual_store_init(dualstore, p);

    /* backpropagate the gradient to grad_pos_attr.[POS], and out.[attribute]. */
    FastPMPainter diffpainter[1];

    /* write the readout result to out*/
    dual_store_set_read_from(dualstore, p);

    /* gradient over the attribute */
    if(attribute > 0) {
        dual_store_set_write_to(dualstore, grad_attr);
        fastpm_readout(painter, y, (FastPMStore*) dualstore, get_position, attribute);
    }

    dual_store_set_write_to(dualstore, grad_pos);
    /* gradient over the particle position. */
    int POS[3] = {PACK_POS_X, PACK_POS_Y, PACK_POS_Z};
    int d;
    for(d = 0; d < 3; d++) {
        fastpm_painter_init_diff(diffpainter, painter, d);
        fastpm_readout(diffpainter, y, (FastPMStore*) dualstore, get_position, POS[d]);
    }
    int i;
    for(i = 0; i < p->np; i ++) {
        double A = attribute?p->to_double(p, i, attribute):1.0;
        for(d = 0; d < 3; d++) {
            grad_pos->x[i][d] *= A;
        }
    }
}

void
fastpm_readout_gradient(FastPMPainter * painter, FastPMStore * y,
    FastPMFloat * canvas, FastPMStore * p, fastpm_posfunc get_position, int attribute,
    FastPMFloat * grad_canvas, FastPMStore * grad_pos)
{
    /* gradient over canvas is stored in out1, gradient over position is stored in out.*/
    DualStore dualstore[1];

    dual_store_init(dualstore, p);

    /* write the readout result to out*/
    /* chain the attribute gradient from y */
    dual_store_set_read_from(dualstore, y);
    dual_store_set_write_to(dualstore, grad_pos);

    /* backpropagate the gradient to out.[POS], and out_canvas */
    FastPMPainter diffpainter[1];

    fastpm_store_copy(p, grad_pos);

    /* gradient over the canvas */
    fastpm_paint(painter, grad_canvas, (FastPMStore*) dualstore, get_position, attribute);

    /* gradient over the particle position. */
    int POS[3] = {PACK_POS_X, PACK_POS_Y, PACK_POS_Z};

    int d;
    for(d = 0; d < 3; d++) {
        fastpm_painter_init_diff(diffpainter, painter, d);
        fastpm_readout(diffpainter, canvas, (FastPMStore*) dualstore, get_position, POS[d]);
    }
    int i;
    for(i = 0; i < p->np; i ++) {
        double v = attribute?y->to_double(y, i, attribute):1.0;
        for(d = 0; d < 3; d++) {
            grad_pos->x[i][d] *= v;
        }
    }

}
