#include <fastpm/libfastpm.h>

#include "pmesh/_window_imp.c"

struct FastPMMeshPainterPrivate {
    PMMeshPainter basepainter;
};


FastPMMeshPainter *
fastpm_mesh_painter_new(
    FastPMMesh * mesh,
    enum FastPMMeshPainterType type,
    int order[])
{
    FastPMMeshPainter * self = malloc(sizeof(self[0]));
    self->type = type;
    self->mesh = mesh;
    int d;
    for(d = 0; d < mesh->ndim; d ++) {
        self->order[d] = order?order[d]:0;
    }
    self->priv = malloc(sizeof(self->priv));
    PMeshPainter painter = {
        .type = _pmesh_painter_type_table[self->type];
        .support = 0,
        .ndim = self->mesh->ndim,
        .canvas_dtype_elsize = sizeof(FastPMFloat),
    };
    int d;
    for(d = 0; d < self->mesh->ndim; d++) {
        painter.order[d] = self->order[d];
        painter.Nmesh[d] = self->mesh->Nmesh[d];
        painter.scale[d] = self->mesh->InvCellSize[d];
        painter.translate[d] = -self->mesh->ral.start[d];
        painter.size[d] = self->ral.shape[d];
        /* convert strides to bytes */
        painter.strides[d] = self->ral.strides[d] * sizeof(FastPMFloat);
    }
    self->priv->basepainter = painter;
    return self;
}

void
fastpm_mesh_painter_free(FastPMMeshPainter * self)
{
    free(self->priv);
    free(self);
}

void
fastpm_mesh_painter_paint(FastPMMeshPainter * self,
    FastPMColumn * x,
    FastPMColumn * m,
    FastPMFloat * canvas
)
{
    PMeshPainter painter = self->priv->basepainter;
    painter.canvas = canvas;
    pmesh_painter_init(&painter);
    double pos[3];
    double mass;
    int i;
    for(i = 0; i < x->size; i ++) {
        fastpm_column_to_double(x, i, pos);
        fastpm_column_to_double(m, i, &mass);
        pmesh_painter_paint(&painter, pos, mass);
    }
}

void
fastpm_mesh_painter_readout(FastPMMeshPainter * self,
    FastPMColumn * x,
    FastPMColumn * m,
    FastPMFloat * canvas
)
{
    PMeshPainter painter = self->priv->basepainter;
    painter.canvas = canvas;
    pmesh_painter_init(&painter);
    double pos[3];
    double mass;
    int i;
    for(i = 0; i < x->size; i ++) {
        fastpm_column_to_double(x, i, pos);
        mass = pmesh_painter_readout(&painter, pos);
        fastpm_column_from_double(m, i, &mass);
    }
}
