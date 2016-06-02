FASTPM_BEGIN_DECLS

typedef enum { FASTPM_PAINTER_CIC, FASTPM_PAINTER_LINEAR, FASTPM_PAINTER_LANCZOS} FastPMPainterType;

struct FastPMPainter {
    PM * pm;

    void   (*paint)(FastPMPainter * painter, FastPMFloat * canvas, double pos[3], double weight);
    double (*readout)(FastPMPainter * painter, FastPMFloat * canvas, double pos[3]);
    fastpm_kernelfunc kernel;

    int support;
    int Npoints; /* (2 * support) ** 3 */
};

void fastpm_painter_init(FastPMPainter * painter, PM * pm,
        FastPMPainterType type, int support);

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

FASTPM_END_DECLS
