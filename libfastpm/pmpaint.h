typedef void   (*fastpm_paintfunc)(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight);
typedef double (*fastpm_readoutfunc)(FastPMPainter * painter, FastPMFloat * canvas, double pos[3]);

struct FastPMPainter {
    PM * pm;

    fastpm_paintfunc paint;
    fastpm_readoutfunc readout;
    fastpm_kernelfunc kernel;

    int support;
    int Npoints;
    int strides[3];
};

void fastpm_painter_init(FastPMPainter * painter, PM * pm,
        fastpm_kernelfunc kernel, int support);

void fastpm_painter_init_cic(FastPMPainter * painter);

void fastpm_painter_paint(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight);
double fastpm_painter_readout(FastPMPainter * painter, FastPMFloat * canvas, double pos[3]);

void
fastpm_paint_store(FastPMPainter * pm, FastPMFloat * canvas,
        PMStore * p, size_t size,
        fastpm_posfunc get_position, int attribute);

void
fastpm_readout_store(FastPMPainter * pm, FastPMFloat * canvas,
        PMStore * p, size_t size,
        fastpm_posfunc get_position, int attribute);

static double
FASTPM_PAINTER_KERNEL_LINEAR(double x, int support) {
    return 1.0 - fabs(x / support);
}
