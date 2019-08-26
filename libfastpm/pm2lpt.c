#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/transfer.h>

#include "pmpfft.h"
#include "pmghosts.h"
#include "pm2lpt.h"

void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, FastPMStore * p, double shift[3], FastPMKernelType type)
{
    /* read out values at locations with an inverted shift */
    int potorder, gradorder, deconvolveorder;
    fastpm_kernel_type_get_orders(type, &potorder, &gradorder, &deconvolveorder);

    /* calculate dx1, dx2, for initial fluctuation delta_k.
     * shift: martin has shift = 0.5, 0.5, 0.5.
     * Use shift of 0, 0, 0 if in doublt. 
     *   */
    ptrdiff_t i;
    int d;

    /* It is important to (de-)shift the particles before creating ghosts 
     * Because we will read out from the (de-)shifted positions.
     * Otherwise the IC will have artifacts along the edges of domains. */
    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] -= shift[d];
        }
    }
    FastPMPainter painter[1];
    fastpm_painter_init(painter, pm, FASTPM_PAINTER_CIC, 0);

    PMGhostData * pgd = pm_ghosts_create(pm, p, p->attributes | COLUMN_DX1 | COLUMN_DX2, painter->support);
    pm_ghosts_send(pgd, COLUMN_POS);

    FastPMFloat * workspace = pm_alloc(pm);
    FastPMFloat * source =  pm_alloc(pm);
    memset(source, 0, sizeof(source[0]) * pm->allocsize);
    memset(workspace, 0, sizeof(workspace[0]) * pm->allocsize);

    FastPMFloat * field[3];

    for(d = 0; d < 3; d++ ) {
        field[d] = pm_alloc(pm);
        memset(field[d], 0, sizeof(field[d][0]) * pm->allocsize);
    }
    FastPMFieldDescr DX1[] = { {COLUMN_DX1, 0}, {COLUMN_DX1, 1}, {COLUMN_DX1, 2}};
    FastPMFieldDescr DX2[] = { {COLUMN_DX2, 0}, {COLUMN_DX2, 1}, {COLUMN_DX2, 2}};
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};

    for(d = 0; d < 3; d++) {

        fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
        fastpm_apply_diff_transfer(pm, workspace, workspace, d);

        pm_c2r(pm, workspace);

        fastpm_readout_local(painter, workspace, p, p->np, DX1[d]);
        fastpm_readout_local(painter, workspace, pgd->p, pgd->p->np, DX1[d]);

    } 

    for(d = 0; d< 3; d++) {
        fastpm_apply_laplace_transfer(pm, delta_k, field[d], potorder);
        fastpm_apply_diff_transfer(pm, field[d], field[d], d);
        fastpm_apply_diff_transfer(pm, field[d], field[d], d);

        pm_c2r(pm, field[d]);
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];

        fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
        fastpm_apply_diff_transfer(pm, workspace, workspace, d1);
        fastpm_apply_diff_transfer(pm, workspace, workspace, d2);

        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] -= workspace[i] * workspace[i];
        }
    } 
    pm_r2c(pm, source, workspace);
    pm_assign(pm, workspace, source);

    for(d = 0; d < 3; d++) {
        /* 
         * We absorb some the negative factor in za transfer to below;
         *
         * */
        fastpm_apply_laplace_transfer(pm, source, workspace, potorder);
        fastpm_apply_diff_transfer(pm, workspace, workspace, d);

        pm_c2r(pm, workspace);

        /* this ensures x = x0 + dx1(t) + dx2(t) */
        fastpm_apply_multiply_transfer(pm, workspace, workspace, 3.0 / 7);

        fastpm_readout_local(painter, workspace, p, p->np, DX2[d]);
        fastpm_readout_local(painter, workspace, pgd->p, pgd->p->np, DX2[d]);
    }

    pm_ghosts_reduce(pgd, COLUMN_DX1, FastPMReduceAddFloat, NULL);
    pm_ghosts_reduce(pgd, COLUMN_DX2, FastPMReduceAddFloat, NULL);

#ifdef PM_2LPT_DUMP
    fwrite(p->dx1, sizeof(p->dx1[0]), p->np, fopen("dx1.f4x3", "w"));
    fwrite(p->dx2, sizeof(p->dx2[0]), p->np, fopen("dx2.f4x3", "w"));
#endif

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += shift[d];
        }
    }

    for(d = 0; d < 3; d ++) {
        pm_free(pm, field[2-d]);
    }
    pm_free(pm, source);
    pm_free(pm, workspace);

    pm_ghosts_free(pgd);

}

// Interpolate position and velocity for snapshot at a=aout
void
pm_2lpt_evolve(double aout, FastPMStore * p, FastPMCosmology * c, int zaonly)
{
    int np = p->np;

    double D1 = GrowthFactor(aout, c);
    double D2 = GrowthFactor2(aout, c);

    double Dv1 = D1 * aout * aout * HubbleEa(aout, c) * DLogGrowthFactor(aout, c);
    double Dv2 = D2 * aout * aout * HubbleEa(aout, c) * DLogGrowthFactor2(aout, c);
    if(zaonly) {
        D2 = 0;
        Dv2 = 0;
    }

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += D1 * p->dx1[i][d] + D2 * p->dx2[i][d];

            if(p->v)
                p->v[i][d] = (p->dx1[i][d]* Dv1 + p->dx2[i][d]*Dv2);
        }
    }
    p->meta.a_x = p->meta.a_v = aout;
}

