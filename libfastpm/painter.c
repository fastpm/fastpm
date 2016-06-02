#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"
#include "pmstore.h"

/* from cic.c */
void fastpm_painter_init_cic(FastPMPainter * painter);

static void
_generic_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight);
static double
_generic_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3]);

static double
_linear_kernel(double x, int support) {
    return 1.0 - fabs(x / support);
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

static double
_lanczos_kernel(double x, int support) {
    if(x >= support || x <= - support) return 0;
    return __sinc__(x) * __sinc__(x / support);
}

void
fastpm_painter_init(FastPMPainter * painter, PM * pm,
    FastPMPainterType type, int support)
{
    painter->pm = pm;
    painter->paint = _generic_paint;
    painter->readout = _generic_readout;
    painter->support = support;

    switch(type) {
        case FASTPM_PAINTER_CIC:
            fastpm_painter_init_cic(painter);
            painter->kernel = NULL;
        break;
        case FASTPM_PAINTER_LINEAR:
            painter->kernel = _linear_kernel;
        break;
        case FASTPM_PAINTER_LANCZOS:
            painter->kernel = _lanczos_kernel;
        break;
    }
    int nmax = 1;
    int d;
    for(d = 0; d < 3; d++) {
        painter->strides[d] = nmax;
        nmax *= (2 * support);
    }
    painter->Npoints = nmax;
}

void
fastpm_painter_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight)
{
    painter->paint(painter, canvas, pos, weight);
}

double
fastpm_painter_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3])
{
    return painter->readout(painter, canvas, pos);
}

static void
_fill_k(FastPMPainter * painter, double pos[3], int ipos[3], float k[3][100])
{
    PM * pm = painter->pm;
    double gpos[3];
    int d;
    for(d = 0; d < 3; d++) {
        gpos[d] = pos[d] * pm->InvCellSize[d];
        ipos[d] = floor(gpos[d]) - (painter->support - 1);
        double dx = gpos[d] - ipos[d];
        int i;
        double sum = 0;
        for(i = 0; i < 2 * painter->support; i ++) {
            k[d][i] = painter->kernel(dx - i, painter->support);
            sum += k[d][i];
        }
        /* normalize the kernel to conserve mass */
        for(i = 0; i < 2 * painter->support; i ++) {
            k[d][i] /= sum;
        }
        ipos[d] -= pm->IRegion.start[d];
    }
}

static void
_generic_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight)
{
    PM * pm = painter->pm;
    int n;
    int d;
    int ipos[3];
    /* the max support is 50 */
    float k[3][100];

    _fill_k(painter, pos, ipos, k);

    for(n = 0; n < painter->Npoints; n ++) {
        float kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n / painter->strides[d]) % (2 * painter->support);

            int targetpos = ipos[d] + rel;
            kernel *= k[d][rel];
            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            ind += pm->IRegion.strides[d] * targetpos;
            if(targetpos < 0 || targetpos >= pm->IRegion.size[d]) 
                goto outside;
        }
#pragma omp atomic
        canvas[ind] += weight * kernel;

    outside:
        continue;
    }
    return;
}

static double
_generic_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3])
{
    PM * pm = painter->pm;
    double value = 0;
    int n;
    int ipos[3];
    float k[3][100];
    int d;

    _fill_k(painter, pos, ipos, k);

    for(n = 0; n < painter->Npoints; n ++) {
        float kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < 3; d++) {
            int rel = (n / painter->strides[d]) % (2 * painter->support);

            kernel *= k[d][rel];

            int targetpos = ipos[d] + rel;

            while(targetpos >= pm->Nmesh[d]) {
                targetpos -= pm->Nmesh[d];
            }
            while(targetpos < 0) {
                targetpos += pm->Nmesh[d];
            }
            if(targetpos >= pm->IRegion.size[d]) {
                goto outside;
            }
            ind += pm->IRegion.strides[d] * targetpos;
        }
        value += kernel * canvas[ind];
outside:
        continue;
    }
    return value;
}

void
fastpm_paint_store(FastPMPainter * painter, FastPMFloat * canvas,
    PMStore * p, size_t size,
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
        fastpm_painter_paint(painter, canvas, pos, weight);
    }
}

void
fastpm_readout_store(FastPMPainter * painter, FastPMFloat * canvas,
    PMStore * p, size_t size,
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
        double weight = fastpm_painter_readout(painter, canvas, pos);
        p->from_double(p, i, attribute, weight);
    }
}
