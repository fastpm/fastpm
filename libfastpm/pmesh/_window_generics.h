static void
mkname(_generic_paint) (PMeshPainter * painter, double pos[], double weight)
{
    int ipos[painter->ndim];
    /* the max support is 32 */
    double k[painter->ndim][64];

    char * canvas = (char*) painter->canvas;

    _fill_k(painter, pos, ipos, k);

    int rel[painter->ndim];
    int d;
    for(d =0; d < painter->ndim; d ++ ) rel[d] = 0;

    int s2 = painter->support;
    while(rel[0] != s2) {
        double kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < painter->ndim; d++) {
            int r = rel[d];
            int targetpos = ipos[d] + r;
            kernel *= k[d][r];
            if(painter->Nmesh[d] > 0) {
                while(targetpos >= painter->Nmesh[d]) {
                    targetpos -= painter->Nmesh[d];
                }
                while(targetpos < 0) {
                    targetpos += painter->Nmesh[d];
                }
            }
            if(UNLIKELY(targetpos >= painter->size[d]))
                goto outside;
            if(UNLIKELY(targetpos < 0))
                goto outside;
            ind += painter->strides[d] * targetpos;
        }
#pragma omp atomic
        * (FLOAT*) (canvas + ind) += weight * kernel;

    outside:
        rel[painter->ndim - 1] ++;
        for(d = painter->ndim - 1; d > 0; d --) {
            if(UNLIKELY(rel[d] == s2)) {
                rel[d - 1] ++;
                rel[d] = 0;
            }
        }
    }
    return;
}

static double
mkname(_generic_readout) (PMeshPainter * painter, double pos[])
{
    double value = 0;
    int ipos[painter->ndim];
    double k[painter->ndim][64];

    char * canvas = (char*) painter->canvas;

    _fill_k(painter, pos, ipos, k);

    int rel[painter->ndim];
    int d;
    for(d =0; d < painter->ndim; d++) rel[d] = 0;

    int s2 = painter->support;
    while(rel[0] != s2) {
        double kernel = 1.0;
        ptrdiff_t ind = 0;
        int d;
        for(d = 0; d < painter->ndim; d++) {
            int r = rel[d];

            kernel *= k[d][r];

            int targetpos = ipos[d] + r;

            if(painter->Nmesh[d] > 0) {
                while(targetpos >= painter->Nmesh[d]) {
                    targetpos -= painter->Nmesh[d];
                }
                while(targetpos < 0) {
                    targetpos += painter->Nmesh[d];
                }
            }
            if(UNLIKELY(targetpos >= painter->size[d])) {
                goto outside;
            }
            if(UNLIKELY(targetpos < 0))
                goto outside;
            ind += painter->strides[d] * targetpos;
        }
        value += kernel * *(FLOAT* )(canvas + ind);
outside:
        rel[painter->ndim - 1] ++;
        for(d = painter->ndim - 1; d > 0; d --) {
            if(UNLIKELY(rel[d] == s2)) {
                rel[d - 1] ++;
                rel[d] = 0;
            }
        }
    }
    return value;
}
